import argparse
import sys
import gzip
import subprocess
import os
import fileinput
import time
from collections import deque
from multiprocessing import Process, Queue

# extend path to include the provided cutadapt
trimPath, trimFile = os.path.split(os.path.realpath(__file__))
cutPath = os.path.join(trimPath, 'miRge.seqUtils', 'cutadapt-1.7.1', 'bin')
env = os.environ.copy()
env['PATH'] = '{0}:{1}'.format(env['PATH'], cutPath)

class GenericIterator(object):
    gz = False
    CHUNK_SIZE = 2**16
    UNCONSUMED = ''
    contents = []

    def __init__(self, filename, **kwrds):
        if isinstance(filename, basestring) and filename.endswith('.gz'):
            self.gz = True
            self.filename = gzip.GzipFile(filename)
        elif isinstance(filename, basestring):
            self.filename = open(filename)
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
        if self.gz:
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
        else:
            return self.filename.next().strip()

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
    def __init__(self, queue=None, cutadapt=None, adapter=None, phred64=False):
        super(Worker, self).__init__()
        self.queue=queue
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
                                  '64' if self.phred64 else '33', '--quiet',
                                  '-o', '{0}.trim{1}'.format(outfile, outext), filename], env=env)
            p.communicate()

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
    chunksize = 1000000
    adapter = args.adapter

    # make the new files, since we don't know its size from the beginning and it'd be wasteful
    # to read it twice, split it to a million reads per file and process as such
    read_queue = Queue()
    workers = []
    for i in xrange(args.threads):
        worker = Worker(queue=read_queue, cutadapt=args.cutadapt, phred64=phred, adapter=adapter)
        workers.append(worker)
        worker.start()

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
                for worker in workers:
                    worker.phred64 = True
        o.write('%s\n' % '\n'.join(reads))

    if index % chunksize:
        o.flush()
        o.close()
        read_queue.put(o.name)
        tmpfiles.append(o.name)

    # poison pill to stop workers
    for i in range(args.threads):
        read_queue.put(None)

    while any([i.is_alive() for i in workers]):
        time.sleep(1)

    # recombine all our files
    with dest as fout:
        for entry in fileinput.input(files):
            fout.write(entry)

    # delete temp files
    for filename in files+tmpfiles:
        os.remove(filename)

    # log
    with open('{0}.log'.format(logfile), 'wb') as o:
        o.write('Starting reads: {0}\n'.format(index+1))
        trimmed = FastqIterator(dest.name)
        for dest_index, dest_read in enumerate(trimmed):
            pass
        o.write('Processed reads: {0}\n'.format(dest_index+1))

    sys.stdout.write('{0}\n'.format(phred))
    return phred

if __name__ == "__main__":
    sys.exit(main())
