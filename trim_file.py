import argparse
import sys
import gzip
import subprocess
import os
import fileinput
import time
import re
from collections import deque
from multiprocessing import Process, Queue
from distutils.spawn import find_executable

# extend path to include the provided cutadapt
trimPath, trimFile = os.path.split(os.path.realpath(__file__))
cutPath = os.path.join(trimPath, 'miRge.seqUtils', 'cutadapt-1.7.1', 'bin')
env = os.environ.copy()
env['PATH'] = '{0}:{1}'.format(env['PATH'], cutPath)

# parse trimmed reads
trimParse = re.compile(r'Trimmed reads\:\s+(?P<trimmed>\d+)')
qtrimParse = re.compile(r'Quality\-trimmed\:\s+(?P<trimmed>\d+)')
tstrimParse = re.compile(r'Too short reads\:\s+(?P<trimmed>\d+)')
processedParse = re.compile(r'Processed reads\:\s+(?P<processed>\d+)')

class GenericIterator(object):
    gz = False
    CHUNK_SIZE = 2**16+8
    UNCONSUMED = ''
    contents = deque()

    def __init__(self, filename, **kwrds):
        if isinstance(filename, basestring) and filename.endswith('.gz'):
            self.gz = True
            name = filename
        elif isinstance(filename, basestring):
            name = filename
        elif isinstance(filename, (file,)):
            if filename.name.endswith('.gz'):
                self.gz = True
            name = filename.name
        else:
            raise TypeError
        if self.gz:
            self.filename = gzip.GzipFile(fileobj=open(name, 'rb', buffering=self.CHUNK_SIZE), mode='rb')
        else:
            self.filename = open(name, mode='rb', buffering=self.CHUNK_SIZE)

    def __iter__(self):
        return self

    def next(self):
        if self.contents:
            return self.contents.popleft()
        new_contents = self.filename.read(self.CHUNK_SIZE)
        if not new_contents:
            if self.UNCONSUMED:
                return self.UNCONSUMED
            raise StopIteration
        if new_contents and new_contents[-1] != '\n':
            new_uc_index = new_contents.rfind('\n')+1
            new_unconsumed = new_contents[new_uc_index:]
            new_contents = new_contents[:new_uc_index]
        else:
            new_unconsumed = ''
        new_contents = self.UNCONSUMED+new_contents
        self.contents = new_contents.split('\n')
        self.contents = filter(None, self.contents)
        self.UNCONSUMED = new_unconsumed
        self.contents = deque(self.contents)
        if self.contents:
            return self.contents.popleft()

class FastqIterator(GenericIterator):
    def __init__(self, filename):
        """
        Optional argument: delimiter -- > default
        """
        super(FastqIterator, self).__init__(filename)

    def next(self):
        # return is sequence header, sequence, quality header, quality sequence
        _next = super(FastqIterator, self).next
        return (_next(), _next(), _next(), _next())


class Worker(Process):
    def __init__(self, queue=None, results=None, cutadapt=None, adapter=None, phred64=False):
        super(Worker, self).__init__()
        self.queue=queue
        self.results = results
        self.phred64 = False
        self.cutadapt = 'cutadapt' if cutadapt is None else cutadapt
        adapter_flag = '-a'
        self.ion = False
        if adapter.startswith('+'):
            adapter_flag = '-u'
            self.ion = True
        if adapter == 'none':
            adapter_command = ['--no-trim']
        else:
            adapter_command = [adapter_flag, adapter]
            if self.ion is False:
                adapter_command.append('--discard-untrimmed')
        self.adapter_command = ' '.join(adapter_command)

    def run(self):
        for filename in iter(self.queue.get, None):
            outfile, outext = os.path.splitext(filename)
            p = subprocess.Popen([self.cutadapt, '-q', '10', '-m', '16', self.adapter_command,
                                  '-e', '0.12', '--quality-base',
                                  '64' if self.phred64 else '33',
                                  '-o', '{0}.trim{1}'.format(outfile, outext), filename], env=env, stdout=subprocess.PIPE)
            sout, serr = p.communicate()
            matched_count = 0
            if self.ion:
                matched = qtrimParse.search(sout)
                if matched:
                    matched_count += int(matched.group('trimmed'))
                matched = tstrimParse.search(sout)
                if matched:
                    matched_count += int(matched.group('trimmed'))
            else:
                matched = trimParse.search(sout)
                if matched:
                    matched_count = matched.group('trimmed')
                else:
                    matched_count = 0
            processed = processedParse.search(sout)
            self.results.put({'trimmed': matched_count, 'processed': processed.group('processed')} if matched_count and processed else None)



parser = argparse.ArgumentParser()
parser.add_argument('--cutadapt', type=str)
parser.add_argument('--adapter', type=str)
parser.add_argument('--infile', type=argparse.FileType('rb'))
parser.add_argument('--outfile', type=argparse.FileType('wb'))
parser.add_argument('--threads', type=int, default=1)
parser.add_argument('--phred64', action='store_true')

def main():
    args = parser.parse_args()
    source = FastqIterator(args.infile)
    dest = args.outfile
    outfile, outext = os.path.splitext(args.infile.name)
    phred = args.phred64 or 33
    threads = args.threads
    gzipped = False
    if outext == '.gz':
        gzipped = True
        outfile, outext = os.path.splitext(outfile)
    logfile = args.infile.name

    # Since we don't know the file sizes from the beginning and it'd be wasteful
    # to read it twice, split it to x reads per file and process as such
    chunksize = 1000000
    adapter = args.adapter

    read_queue = Queue()
    result_queue = Queue()

    workers = []
    for i in xrange(threads):
        worker = Worker(queue=read_queue, results=result_queue, cutadapt=args.cutadapt, phred64=phred==64, adapter=adapter)
        workers.append(worker)
        worker.start()

    o = None
    files = []
    tmpfiles = []
    if threads == 1:
        filename = args.infile.name
        read_queue.put(filename)
        fbase, fext = os.path.splitext(filename)
        files.append('{0}.trim{1}'.format(fbase, fext))
    else:
        for index, reads in enumerate(source):
            newfile = index % chunksize == 0
            if newfile:
                if o is not None:
                    o.flush()
                    o.close()
                    read_queue.put(o.name)
                    tmpfiles.append(o.name)
                filename = '{0}_{1}{2}'.format(outfile, str(index/chunksize), outext)
                fbase, fext = os.path.splitext(filename)
                files.append('{0}.trim{1}'.format(fbase, fext))
                o = open(filename, 'wb')
            if index < 1000 and phred == 33:
                if any([i for i in reads[3] if ord(i) > 74]):
                    phred = 64
            o.write('%s\n' % '\n'.join(reads))

        if index % chunksize:
            o.flush()
            o.close()
            read_queue.put(o.name)
            tmpfiles.append(o.name)

    # poison pill to stop workers
    for i in xrange(threads):
        read_queue.put(None)

    while any([i.is_alive() for i in workers]):
        time.sleep(0.01)

    # recombine all our files
    if threads != 1:
        with dest as fout:
            # see if we can use cat, else use python
            cat = find_executable('cat')
            if cat is not None:
                cmd = [cat]
                cmd += files
                p = subprocess.Popen(cmd, stdout=fout, env=env)
                p.communicate()
            else:
                for entry in fileinput.input(files):
                    fout.write(entry)


        # delete temp files
        for filename in files+tmpfiles:
            os.remove(filename)

    # log, if the result queue has the trimmed count use it, else figure it out. We do it this way incase cutadapt changes
    # their output
    # poison pill for results
    result_queue.put(-1)

    with open('{0}.log'.format(logfile), 'wb') as o:
        results = [result for result in iter(result_queue.get, -1)]
        if None in results:
            # something changed with cutadapt
            if threads == 1:
                for index, reads in enumerate(source):
                    pass
            o.write('Starting reads: {0}\n'.format(index+1))
            trimmed = FastqIterator(dest.name)
            for dest_index, dest_read in enumerate(trimmed):
                pass
            o.write('Processed reads: {0}\n'.format(dest_index+1))
        else:
            o.write('Starting reads: {0}\n'.format(sum(map(int,[i.get('processed',0) for i in results]))))
            o.write('Processed reads: {0}\n'.format(sum(map(int,[i.get('trimmed',0) for i in results]))))


    sys.stdout.write('{0}\n'.format(phred))
    return phred

if __name__ == "__main__":
    sys.exit(main())
