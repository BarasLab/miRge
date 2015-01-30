import argparse, sys, gzip, subprocess, os, fileinput, time
from collections import deque
from multiprocessing import Process, Queue

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
    def __init__(self, queue=None, cutadapt=None, adapter='TGGAATTCTCGGGTGCCAAGGAACTCCAG'):
        super(Worker, self).__init__()
        self.queue=queue
        self.cutadapt = cutadapt
        self.adapter = adapter

    def run(self):
        for filename in iter(self.queue.get, None):
            p = subprocess.Popen([self.cutadapt, '-q', '10', '-m', '16', '-a',
                                  self.adapter, '-e', '0.12', '--quiet',
                                  '--discard-untrimmed',
                                  '-o', '{0}_trim'.format(filename), filename])
            p.communicate()

parser = argparse.ArgumentParser()
parser.add_argument('--cutadapt', type=str)
parser.add_argument('--infile', type=argparse.FileType('rb'))
parser.add_argument('--threads', type=int, default=1)

def main():
    args = parser.parse_args()
    source = FastqIterator(args.infile)
    outfile, outext = os.path.splitext(args.infile.name)
    if outext == '.gz':
        outfile, outext = os.path.splitext(outfile)
    chunksize = 1000000
    # make the new files, since we don't know its size from the beginning and it'd be wasteful
    # to read it twice, split it to a million reads per file and process as such
    read_queue = Queue()
    workers = []
    for i in xrange(args.threads):
        worker = Worker(queue=read_queue, cutadapt=args.cutadapt)
        workers.append(worker)
        worker.start()

    o = None
    files = []
    for index, reads in enumerate(source):
        newfile = index % chunksize == 0
        if newfile:
            if o is not None:
                o.flush()
                o.close()
                read_queue.put(o.name)
            filename = '{0}_{1}{2}'.format(outfile, str(index/chunksize), outext)
            files.append('{0}_trim'.format(filename))
            o = open(filename, 'wb')
        o.write('%s\n' % '\n'.join(reads))

    if index % chunksize:
        o.flush()
        o.close()
        read_queue.put(o.name)

    # poison pill to stop workers
    for i in range(args.threads):
        read_queue.put(None)

    while any([i.is_alive() for i in workers]):
        time.sleep(1)

    # recombine all our files
    with open('{0}_trim{1}'.format(outfile, outext), 'wb') as fout:
        for entry in fileinput.input(files):
            fout.write(entry)

if __name__ == "__main__":
    sys.exit(main())