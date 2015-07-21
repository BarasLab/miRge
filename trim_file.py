import argparse
import sys
import gzip
import subprocess
import os
import fileinput
import time
from collections import deque
from multiprocessing import Process, Queue
from distutils.spawn import find_executable
from cutadapt.adapters import Adapter, gather_adapters
from cutadapt.scripts.cutadapt import AdapterCutter
from cutadapt.modifiers import QualityTrimmer, UnconditionalCutter
from cutadapt.seqio import Sequence, FastqReader
from cStringIO import StringIO

# extend path to include the provided cutadapt
trimPath, trimFile = os.path.split(os.path.realpath(__file__))
cutPath = os.path.join(trimPath, 'miRge.seqUtils', 'cutadapt-1.7.1', 'bin')
env = os.environ.copy()
env['PATH'] = '{0}:{1}'.format(env['PATH'], cutPath)

class Worker(Process):
    def __init__(self, queue=None, results=None, adapter=None, phred64=False):
        super(Worker, self).__init__()
        self.queue=queue
        self.results = results
        self.phred = 64 if phred64 else 33
        self.modifiers = [QualityTrimmer(0, 10, self.phred)]
        self.adapters = []
        self.error_rate = 0.12
        self.min_length = 16
        if adapter.startswith('+'):
            self.modifiers.append(UnconditionalCutter(int(adapter)))
        elif adapter == 'none':
            self.adapter = None
        else:
            name, seq, where = gather_adapters([adapter], [], []).next()
            self.adapters = [Adapter(seq, where, self.error_rate, name=name)]
            adapter_cutter = AdapterCutter(self.adapters)
            self.modifiers.append(adapter_cutter)

    def run(self):
        # we can't use the sentinel iter(self.queue.get, None) because of some issue with Sequence classes
        sequence = self.queue.get()
        while sequence is not None:
            read = sequence
            sequence = self.queue.get()
            for modifier in self.modifiers:
                read = modifier(read)
            if len(read.sequence) < self.min_length:
                continue
            self.results.put(read)

class Writer(Process):
    def __init__(self, queue=None, outfile=None):
        super(Writer, self).__init__()
        self.queue = queue
        self.outfile = outfile
        self.kept = 0
        # we batch our writes to minimize disk seek
        self.to_write = []

    def run(self):
        read = self.queue.get()
        while read is not None:
            self.to_write.append(read)
            if len(self.to_write) > 100000:
                [i.write(self.outfile) for i in self.to_write]
                self.to_write = []
            self.kept += 1
            read = self.queue.get()
        [i.write(self.outfile) for i in self.to_write]
        self.queue.put(self.kept)

parser = argparse.ArgumentParser()
parser.add_argument('--cutadapt', type=str)
parser.add_argument('--adapter', type=str)
parser.add_argument('--infile', type=argparse.FileType('rb'))
parser.add_argument('--outfile', type=argparse.FileType('wb'))
parser.add_argument('--threads', type=int, default=1)
parser.add_argument('--phred64', action='store_true')

def main():
    args = parser.parse_args()
    dest = args.outfile
    phred = args.phred64 or 33
    threads = args.threads
    logfile = args.infile.name

    adapter = args.adapter

    read_queue = Queue()
    result_queue = Queue()

    workers = []
    for i in xrange(threads):
        worker = Worker(queue=read_queue, results=result_queue, phred64=phred==64, adapter=adapter)
        workers.append(worker)
        worker.start()

    writer = Writer(queue=result_queue, outfile=dest)
    writer.start()

    for index, read in enumerate(FastqReader(args.infile)):
        read_queue.put(read)
    processed = index+1

    # poison pill to stop workers
    for i in xrange(threads):
        read_queue.put(None)

    while any([i.is_alive() for i in workers]):
        time.sleep(0.01)

    # poison pill for results
    result_queue.put(None)

    while writer.is_alive():
        time.sleep(0.01)

    kept_reads = result_queue.get()

    with open('{0}.log'.format(logfile), 'wb') as o:
        o.write('Starting reads: {0}\n'.format(processed))
        o.write('Processed reads: {0}\n'.format(kept_reads))

    sys.stdout.write('{0}\n'.format(phred))
    return phred

if __name__ == "__main__":
    sys.exit(main())
