import os
import subprocess as sup


########################################
######                            ######
######      HELPER FUNCTIONS      ######
######                            ######
########################################


######      HELPER FUNCTION FOR SPLIT-GENOME

def assert_directory(path):
    """
        Asserts directory exists
    """
    path = path.split('/')
    file_exists = os.path.exists('/'.join(path))
    if file_exists:
        print('File exists')
    else:
        print('File does not exists!!!!')
        print('Creating directory and / or file')
        if len(path) > 1:
            directory = ''.join(path[:-1])
            if not os.path.exists(directory):
                os.makedirs(directory)
            else:
                pass
            path = '/'.join([directory, path[-1]])
            with open(path, 'w') as fh:
                fh.close()
        else:
            with open(''.join(path), 'w') as fh:
                fh.close()

def __load_seqs(filename):
    """

    """
    FHIN = open(filename, "r")

    line = FHIN.readline()
    seqs = {}
    while line:
        line = line.strip("\n")
        if line[0] == ">":
            header = line
            seqs[header] = []
        else:
            seqs[header].append(line)

        line = FHIN.readline()
    out_seqs = {}
    for head, seq in seqs.items():
        out_seqs[head] = "".join(seq)

    FHIN.close()
    seqs = None
    return out_seqs


def load_lines(infile):
    """
        Loads lines from file and stores them in a list
    """
    lines = []
    with open(infile, 'r') as fh:
        for line in fh:
            lines.append(line.strip())
        
        fh.close()

    return lines

def __split(seqs, size, step):
    """

    """
    out_seqs = {}
    seq_counter = 0
    for head, seq in seqs.items():
        lseq = len(seq)
        i = 0
        start = i
        while start < lseq:
            end = start + size - 1
            out_seqs[">fragment_" + str(seq_counter)] = seq[start:end]
            start = start + step
            i += 1
            seq_counter += 1

    return out_seqs

def __write_seqs(seqs, outfile):
    """

    """
    FHOUT = open(outfile, "w+")
    
    big_string = ""
    for head, seq in seqs.items():
        big_string += head + "\n" + seq + "\n"

    FHOUT.write(big_string)
    FHOUT.close()

##################################################

######      HELPER FUNCTION FOR SELECT-SEQUENCES

def __load_headers(hd_file):
    """
        Load the headers of the sequences that hitted a subject in genomes database.
    """
    
    FHIN = open(hd_file, "r")

    line = FHIN.readline()
    heads = []
    while line:

        line = line.strip("\n")
        heads.append(line)
        line = FHIN.readline()

    FHIN.close()
    return heads


def __select_seqs(seqs, heads):
    """
        Selects those sequences which header is not in heads list.
    """

    lsq = len(seqs)
    print(str(lsq) + ' current sequences')
    # Sort heads
    heads = list(set(heads))
    lhd = len(heads)
    print(str(lhd) + ' dropped secuences')
    if lhd <= 0:
        return seqs
    list_seqs_keys = list(seqs.keys())
    # Sort list_seqs_keys
#    print(len(list_seqs_keys))
    out_seqs = {}

    for x in list_seqs_keys:
        if x[1:] not in heads:
            out_seqs[x] = seqs[x]

    print(str(len(out_seqs)) + ' maintained sequences')
    return out_seqs

#####################################################################

def check_path_exists(path):
    """
        Checks path exist
    """
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        print(path + ' folder exists!')


def copy_file(source, destiny):
    """
        Copies a file
    """
    p = sup.Popen(['cp', source, destiny], stdout = sup.PIPE, stderr = sup.PIPE)
    p.communicate()
    p.kill()
    if os.path.exists(destiny):
        print(source + ' was copied to ' + destiny)
    else:
        print('Copy process | Error!!!')

def rm_file(path):
    """
        Removes the file specifyied in path
    """
    if os.path.exists(path):
        os.remove(path)
    else:
        print(path + ' not found.')


######      HELPER FUNCTIONS FOR MAPPPING


def __ht2idx_ready():
    """
        Checks if Hisat2 database is allready formated.
    """

    for fname in os.listdir('ht2_idx'):
        if fname.endswith('.1.ht2'):
            return True
        
    return False


def __index(ref_genome, fname_suffix):
    """
        Prepare an index of reference genome (a.k.a. query genome)
    """

    p = sup.Popen(['hisat2-build', \
                ref_genome, \
                'ht2_idx/' + fname_suffix], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)
    p.communicate()
    p.kill()


def __map(ref_genome, input_split, fname_suffix):
    """
        Maps sequences to reference genome
    """

    pid = input_split.split('.')[-2]
    mapping_out_sam = 'mapping/nal_frags.' + pid + '.sam'
    mapping_out_bam = 'mapping/nal_frags.' + pid + '.bam'
    mapping_out_sbam = 'mapping/nal_frags.' + pid + '.sort.bam'
    p = sup.Popen(['hisat2', \
                '-x', 'ht2_idx/' + fname_suffix, \
                '--no-temp-splicesite', \
                '--no-spliced-alignment', \
                '-f', \
                '-U', input_split, \
                '-S', mapping_out_sam], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)

    p.communicate()
    p.kill()
    p = sup.Popen(['samtools', 'view', \
                '-b', '-h', \
                '-f', '0x0', '-F', '0x100', \
                '-o', mapping_out_bam, \
                mapping_out_sam], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)
    p.communicate()
    p.kill()
    p = sup.Popen(['samtools', 'sort', \
                '-o', mapping_out_sbam, \
                '-O', 'BAM', \
                '--reference', ref_genome, \
                mapping_out_bam], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)
    p.communicate()
    p.kill()

#########################################################

######      HELPER FUNCTION FOR ASSEMBLY

def __do_assembly(pid):
    p = sup.Popen(['stringtie', 'mapping/nal_frags.' + str(pid) + '.sort.bam', \
                    '-L', '-s', '1.5', '-g', '50', \
                    '-o', 'gtfs/nal_frags.' + str(pid) + '.gtf'], \
                    stdout = sup.PIPE, \
                    stderr = sup.PIPE)

    p.communicate()
    p.kill()


##########################################################


######################################
######                          ######
######      MAIN FUNCTIONS      ######
######                          ######
######################################


def split_genome(in_file, size, step_size, out_file):
    """
        Split query genome in fragments of determined size and at each determined step
    """
    print("Loading genome...")
    seqs = __load_seqs(in_file)
    print("Splitting genome...")
    seqs = __split(seqs, size, step_size)
    print("Writting to file...")
    __write_seqs(seqs, out_file)

    print(f"Spliting genome {in_file}")


def do_blast(query, db_file, out_blast, evalue, idt, qcov, task):
    """
        Performs a BLASTn task.
    """
    
    p = sup.Popen(['blastn', \
                    '-query', query, \
                    '-db', db_file, \
                    '-out', out_blast, \
                    '-outfmt', '6 qseqid', \
                    '-task', task, \
                    '-perc_identity', idt, \
                    '-qcov_hsp_perc', qcov, \
                    '-evalue', evalue, \
                    '-max_target_seqs', str(1)], \
                    stdout = sup.PIPE, \
                    stderr = sup.PIPE)

    p.communicate()
    p.kill()

def select_sequences(in_file, hd_file):
    """
        Selects those sequences that hit a subject in blast searching
    """
    
    seqs = __load_seqs(in_file)
    
    heads = __load_headers(hd_file)

    seqs = __select_seqs(seqs, heads)

    __write_seqs(seqs, in_file)


def mapping(ref_genome, input_split):
    """

    """
    fname_suffix = ref_genome.split('/')[-1]
    fname_suffix = '.'.join(fname_suffix.split('.')[:-1])

    if not __ht2idx_ready():
        print('ht2 index not found.')
        print('Indexing.')
        __index(ref_genome, fname_suffix)
    else:
        print('ht2 index found.')
    
    print('Mapping to reference genome.')
    __map(ref_genome, input_split, fname_suffix)
    print('Mapping finished.')



def assembly(pid):
    """
        Assembles fragments of query genome using the genome-guided procedure.
    """

    __do_assembly(pid)

