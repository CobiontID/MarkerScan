import argparse
import os 

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action='store', dest='fasta', metavar='FASTA',help='define fasta file')
parser.add_argument("-o", type=str, action='store', dest='outdir', metavar='OUTDIR',help='define outdir')
parser.add_argument("-s", type=int, action='store', dest='seqs', metavar='SEQS',help='define number of seqs')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

#from https://github.com/peterjc/pico_galaxy/blob/e8a3b88c1d33c75b99edc17c27b9fc4d1260696e/tools/protein_analysis/seq_analysis_utils.py

def fasta_iterator(filename, max_len=None, truncate=None):
    """Simple FASTA parser yielding tuples of (title, sequence) strings."""
    handle = open(filename)
    title, seq = "", ""
    for line in handle:
        if line.startswith(">"):
            if title:
                if truncate:
                    seq = seq[:truncate]
                if max_len and len(seq) > max_len:
                    raise ValueError("Sequence %s is length %i, max length %i"
                                     % (title.split()[0], len(seq), max_len))
                yield title, seq
            title = line[1:].rstrip()
            seq = ""
        elif title:
            seq += line.strip()
        elif not line.strip() or line.startswith("#"):
            # Ignore blank lines, and any comment lines
            # between records (starting with hash).
            pass
        else:
            handle.close()
            raise ValueError("Bad FASTA line %r" % line)
    handle.close()
    if title:
        if truncate:
            seq = seq[:truncate]
        if max_len and len(seq) > max_len:
            raise ValueError("Sequence %s is length %i, max length %i"
                             % (title.split()[0], len(seq), max_len))
        yield title, seq
    return


def split_fasta(input_filename, output_filename_base, n, truncate=None, keep_descr=False, max_len=None):
    """Split FASTA file into sub-files each of at most n sequences.
    Returns a list of the filenames used (based on the input filename).
    Each sequence can also be truncated (since we only need the start for
    SignalP), and have its description discarded (since we don't usually
    care about it and some tools don't like very long title lines).
    If a max_len is given and any sequence exceeds it no temp files are
    created and an exception is raised.
    """
    iterator = fasta_iterator(input_filename, max_len, truncate)
    files = []
    try:
        while True:
            records = []
            for i in range(n):
                try:
                    records.append(next(iterator))
                except:
                    break
            if not records:
                break
            new_filename = "%s.%i.fa" % (output_filename_base, len(files))
            handle = open(new_filename, "w")
            if keep_descr:
                for title, seq in records:
                    handle.write(">%s\n" % title)
                    for i in range(0, len(seq), 60):
                        handle.write(seq[i:i + 60] + "\n")
            else:
                for title, seq in records:
                    handle.write(">%s\n" % title.split()[0])
                    for i in range(0, len(seq), 60):
                        handle.write(seq[i:i + 60] + "\n")
            handle.close()
            files.append(new_filename)
            # print "%i records in %s" % (len(records), new_filename)
    except ValueError:
        # Max length failure from parser - clean up
        try:
            handle.close()
        except Exception:
            pass
        for f in files:
            if os.path.isfile(f):
                os.remove(f)
        raise err
    for f in files:
        assert os.path.isfile(f), "Missing split file %r (!??)" % f
    return files

filenames=split_fasta(results.fasta,results.outdir+'/kraken.tax',results.seqs,truncate=None, keep_descr=False, max_len=None)         