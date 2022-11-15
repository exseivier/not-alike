#!/usr/bin/env python3
import random as RND
import click
import utils as CMD


@click.group()
def main():
    """
        MAIN FUNCTION
    """
    pass

@main.command()
@click.option('-g', '--genome', help = 'Query genome FASTA file name', required = True, type = str)
#@click.option('-sof', '--split-outfile', help = 'Query split-genome output FASTA file name', required = True, type = str)
@click.option('-ws', '--window-size', help = 'Fragments from split genome will size equal to window size', required = True, type = int)
@click.option('-ss', '--step-size', help = 'It is the step size in nucleotides the window takes between each cut along the genome', required = True, type = int)
@click.option('-db', '--database-file', help = 'A file name of a file that contains the pathway to BLAST_DB ver4 (*.db) files separated by return character', required = True, type = str)
@click.option('-e', '--evalue', help = 'E-value cutoff', required = True, type = str)
@click.option('-i', '--identity', help = 'Identity percentage cutoff', required = True, type = str)
@click.option('-q', '--qcov', help = 'HSP query coverage cutoff', required = True, type = str)
@click.option('-t', '--task', help = 'BLAST task [blastn | megablast | dc-megablast]', required = True, type = str)
def search(genome, window_size, step_size, database_file, evalue, identity, qcov, task):
    """
        Searches for not alike fragments in query genome
    """
    PID = RND.randrange(1, 9999999999)
    out_split = '.'.join(genome.split('.')[:-1]) + '_split_' + str(window_size) + '_' + str(step_size) + '.fasta'
    out_split = ''.join(out_split.split('/')[-1])
    input_split = 'input_split.' + str(PID) + '.fasta'
    print(out_split + ' was loaded.')
    print(input_split + ' was loaded.')

    CMD.check_path_exists('split_out')
    CMD.check_path_exists('blast_db')
    CMD.check_path_exists('blast_out')
    CMD.check_path_exists('ht2_idx')
    CMD.check_path_exists('mapping')
    CMD.check_path_exists('gtfs')

    out_split = 'split_out/' + out_split
    input_split = 'split_out/' + input_split
    if not CMD.os.path.exists(out_split):
        CMD.split_genome(genome, window_size, step_size, out_split)
    else:
        print('A split-genome file was found!!!')

    CMD.copy_file(out_split, input_split)

    db_files = CMD.load_lines(database_file)

    for f in db_files:
        dbf = '.'.join(f.split('.')[:-1]) + '.db'
        print('Blasting ' + dbf + ' ...')
        CMD.do_blast(input_split, 'blast_db/' + dbf, 'blast_out/out.blast', evalue, identity, qcov, task)
        print('Updating input_split')
        CMD.select_sequences(input_split, 'blast_out/out.blast')

    print('BLASTn searching done!')
    print('Mapping on process.')
    
    CMD.mapping(genome, input_split)

    print('Assembly on process.')
    CMD.assembly(PID)

    print('not-alike has finished.')


if __name__ == '__main__':
    main()
