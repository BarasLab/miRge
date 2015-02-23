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

# extend path to include the provided cutadapt
trimPath, trimFile = os.path.split(os.path.realpath(__file__))
cutPath = os.path.join(trimPath, 'miRge.seqUtils', 'cutadapt-1.7.1', 'bin')
env = os.environ.copy()
env['PATH'] = '{0}:{1}'.format(env['PATH'], cutPath)

# parse trimmed reads
trimParse = re.compile(r'Trimmed reads\:\s+(?P<trimmed>\d+)')

class GenericIterator(object):
    gz = False
    CHUNK_SIZE = 2**16+8
    UNCONSUMED = ''
    contents = []

    def __init__(self, filename, **kwrds):
        if isinstance(filename, basestring) and filename.endswith('.gz'):
            self.gz = True
            self.filename = gzip.GzipFile(filename)
        elif isinstance(filename, basestring):
            self.filename = open(filename, 'rb', buffering=self.CHUNK_SIZE)
        elif isinstance(filename, (file,)):
            if filename.name.endswith('.gz'):
                self.gz = True
                self.filename = gzip.GzipFile(filename.name)
            else:
                self.filename = filename
        else:
            raise TypeError

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
        self.cutadapt = cutadapt
        adapter_flag = '-a'
        if adapter.startswith('+'):
            adapter_flag = '-u'
        if adapter == 'none':
            adapter_command = ['--no-trim']
        else:
            adapter_command = [adapter_flag, adapter, '--discard-untrimmed']
        self.adapter_command = ' '.join(adapter_command)

    def run(self):
        for filename in iter(self.queue.get, None):
            outfile, outext = os.path.splitext(filename)
            p = subprocess.Popen([self.cutadapt, '-q', '10', '-m', '16', self.adapter_command,
                                  '-e', '0.12', '--quality-base',
                                  '64' if self.phred64 else '33',
                                  '-o', '{0}.trim{1}'.format(outfile, outext), filename], env=env, stdout=subprocess.PIPE)
            sout, serr = p.communicate()
            matched = trimParse.search(sout)
            self.results.put(matched.group('trimmed') if matched else None)



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
    gzipped = False
    if outext == '.gz':
        gzipped = True
        logfile = outfile
        outfile, outext = os.path.splitext(outfile)
    else:
        logfile = args.infile.name

    # Since we don't know the file sizes from the beginning and it'd be wasteful
    # to read it twice, split it to x reads per file and process as such
    chunksize = 200000
    adapter = args.adapter

    read_queue = Queue()
    result_queue = Queue()

    o = None
    open_func = gzip.open if gzipped else open
    files = []
    tmpfiles = []
    for index, reads in enumerate(source):
        newfile = index % chunksize == 0
        if newfile:
            if o is not None:
                o.flush()
                o.close()
                read_queue.put(o.name)
                tmpfiles.append(o.name)
            filename = '{0}_{1}{2}{3}'.format(outfile, str(index/chunksize), outext, '.gz' if gzipped else '')
            fbase, fext = os.path.splitext(filename)
            files.append('{0}.trim{1}'.format(fbase, fext))
            o = open_func(filename, 'wb')
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
    for i in range(args.threads):
        read_queue.put(None)
    # poison pill for results
    result_queue.put(-1)

    workers = []
    for i in xrange(args.threads):
        worker = Worker(queue=read_queue, results=result_queue, cutadapt=args.cutadapt, phred64=phred==64, adapter=adapter)
        workers.append(worker)
        worker.start()

    while any([i.is_alive() for i in workers]):
        time.sleep(1)

    # recombine all our files
    with dest as fout:
        for entry in fileinput.input(files):
            fout.write(entry)

    # delete temp files
    for filename in files+tmpfiles:
        os.remove(filename)

    # log, if the result queue has the trimmed count use it, else figure it out. We do it this way incase cutadapt changes
    # their output

    with open('{0}.log'.format(logfile), 'wb') as o:
        results = [result for result in iter(result_queue.get, -1)]
        o.write('Starting reads: {0}\n'.format(index+1))
        if None in results:
            # something changed with cutadapt
            trimmed = FastqIterator(dest.name)
            for dest_index, dest_read in enumerate(trimmed):
                pass
            o.write('Processed reads: {0}\n'.format(dest_index+1))
        else:
            o.write('Processed reads: {0}\n'.format(sum(results)))


    sys.stdout.write('{0}\n'.format(phred))
    return phred

if __name__ == "__main__":
    sys.exit(main())
