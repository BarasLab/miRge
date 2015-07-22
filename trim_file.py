import argparse
import sys
from multiprocessing import Process, Queue
from cutadapt.adapters import Adapter, gather_adapters
from cutadapt.scripts.cutadapt import AdapterCutter
from cutadapt.modifiers import QualityTrimmer, UnconditionalCutter
from cutadapt.seqio import Sequence, FastqReader
from cStringIO import StringIO
#from profilestats import profile

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

    #@profile(print_stats=20, dump_stats=True)
    def run(self):
        # we can't use the sentinel iter(self.queue.get, None) because of some issue with Sequence classes
        results = self.results
        get_func = self.queue.get
        reads = get_func()
        modifiers = self.modifiers
        min_length = self.min_length
        while reads is not None:
            result_batch = []
            for read in reads:
                for modifier in modifiers:
                    read = modifier(read)
                if len(read.sequence) >= min_length:
                    io = StringIO()
                    read.write(io)
                    io.seek(0)
                    result_batch.append(io.read())
                    io.close()
            results.put(result_batch)
            reads = get_func()

class Writer(Process):
    def __init__(self, queue=None, trimmed=None, outfile=None):
        super(Writer, self).__init__()
        self.queue = queue
        self.trimmed = trimmed
        self.outfile = outfile

    #@profile(print_stats=20, dump_stats=True)
    def run(self):
        get_func = self.queue.get
        reads = get_func()
        outfile = self.outfile
        kept = 0
        while reads is not None:
            for read in reads:
                outfile.write(read)#read.write(outfile)
                kept += 1
            reads = get_func()
        self.trimmed.put(kept)

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
    trimmed_queue = Queue()

    workers = []
    for i in xrange(threads):
        worker = Worker(queue=read_queue, results=result_queue, phred64=phred==64, adapter=adapter)
        workers.append(worker)
        worker.start()

    writer = Writer(queue=result_queue, trimmed=trimmed_queue, outfile=dest)
    writer.start()

    batch = []
    for index, read in enumerate(FastqReader(args.infile)):
        batch.append(read)
        if index % 10000 == 0:
            read_queue.put(batch)
            batch = []
    read_queue.put(batch)
    processed = index+1

    # poison pill to stop workers
    for i in xrange(threads):
        read_queue.put(None)

    for i in workers:
        i.join()

    # poison pill for writers
    result_queue.put(None)

    # wait for writing to finish
    writer.join()

    trimmed_queue.put(None)

    dest.flush()
    dest.close()

    kept_reads = sum([i for i in iter(trimmed_queue.get, None)])

    with open('{0}.log'.format(logfile), 'wb') as o:
        o.write('Starting reads: {0}\n'.format(processed))
        o.write('Processed reads: {0}\n'.format(kept_reads))

    sys.stdout.write('{0}\n'.format(phred))
    return phred

if __name__ == "__main__":
    sys.exit(main())
