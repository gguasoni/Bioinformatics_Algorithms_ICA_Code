from Bio import SeqIO
import argparse, os
import programme_settings as ps

################################################
# Part 1: Setting up the Commandline Interface #
################################################

ps.read()

def argument_taker():
    parser = argparse.ArgumentParser(prog = "BA_ICA", formatter_class=argparse.RawDescriptionHelpFormatter, description = """\
             Welcome to the BLAST101 Search Interface
             Analyze your sequence using:
                - BLAST101
                - Smith-Waterman
                - Statistical analysis""")

    #Database Input
    parser.add_argument('--database', required = True, help = "Enter database sequence file")

    #Query Sequence Input
    parser.add_argument('--sequence', required = True, help = "The input sequence")

    #Allowing the user to pick between algorithms and statistical analysis
    parser.add_argument('--method', default = 'blast101', required = True, choices=['blast101', 'smith-waterman', 'statistical_analysis'], help = "If you would like to run blast or smith-waterman, select the relevant option. If you would like to conduct statistical analysis, select that option.")

    #Allowing the user to run the tests to check that the UI works
    parser.add_argument('--tests', default = 'false', choices=['true', 'false'], help = "Using this argument will allow you to run tests that check if this commandline works.")

    #abbreviating to args
    args = parser.parse_args()

    database_path = args.database
    sequence_path = args.sequence
    algo = args.method
    run_tests2 = args.tests



    question = input(f"You have entered {database_path} as your database sequence and {sequence_path} as your query sequence, using {algo} as your preferred method of analysis.\n Do you wish to proceed with these files? y/n")
    if question.lower().startswith("y"):
        print("Proceeding with sequences.")

    elif question.lower().startswith("n"):
        raise ValueError("Please input your desired sequences.")

    else:
        raise ValueError("Response unintelligible.")

    # Checking that the files exist and are populated
    does_it_exist(database_path, sequence_path)
    #Checking that the files are FASTAs
    but_does_it_fasta(database_path, sequence_path)
    #Checking that the proteins
    protein_or_not(database_path, sequence_path)

    sequencer = next(SeqIO.parse(args.sequence, "fasta"))
    args.sequence = str(sequencer.seq)

    return args


#Checking that the referenced sequences exist
def does_it_exist(database_path, sequence_path):
    if os.path.exists(database_path) and os.path.exists(sequence_path):
        if os.path.getsize(database_path) and os.path.getsize(sequence_path) != 0:
            print(f"{database_path} and {sequence_path} exist and are populated. Proceeding with file format validation...")
            return True
        else:
            raise ValueError(f"Both {database_path} and {sequence_path} exist but at least one is empty. Please check your files.")
    else:
        raise FileNotFoundError(f"Could not find at least one of {database_path} and {sequence_path}. Please check your directory.")


def but_does_it_fasta(database_path, sequence_path):
    if database_path.lower().endswith((".fasta", ".fa")) and sequence_path.lower().endswith((".fasta", ".fa")):
        print("Input is in the correct format. Proceeding with input validation...")
        return True
    else:
        raise ValueError("This system only accepts FASTA files. Please check your input files.")

#Validating that the sequences are protein sequences
def protein_or_not(database_path, sequence_path):
    prot_alph = set("ACDEFGHIKLMNPQRSTVWY")
    dna_alph = set("ACGT")

    #checking the database sequence first
    database_check = next(SeqIO.parse(database_path, "fasta"))
    database_check = str(database_check.seq).upper()

    db_protein_count = 0
    db_dna_count = 0

    for char in database_check:
        if char not in prot_alph:
            raise ValueError(f"This sequence contains {char}, which is not an amino acid.\nThis system only accepts protein sequences.")
        if char in prot_alph:
            db_protein_count += 1
        if char in dna_alph:
            db_dna_count += 1

    sequence_check = next(SeqIO.parse(sequence_path, "fasta"))
    sequence_check = str(sequence_check.seq).upper()

    sequence_protein_count = 0
    sequence_dna_count = 0

    for char in sequence_check:
        if char not in prot_alph:
            raise ValueError(f"This sequence contains {char}, which is not an amino acid.\nThis system only accepts protein sequences.")

        if char in prot_alph:
            sequence_protein_count += 1
        if char in dna_alph:
            sequence_dna_count += 1


    if db_protein_count > db_dna_count and sequence_protein_count > sequence_dna_count:
        print(f"{database_path} and {sequence_path} are protein sequences.\nSequences accepted.")
        return True

    elif (db_protein_count == db_dna_count) or (sequence_protein_count == sequence_dna_count):
        query = input("At least one of these sequences only contains amino acid residues with the same codes as DNA nucleotides (i.e., 'A', 'C', 'G', 'T'). Are you sure they are both protein sequences? y/n")
        if query.lower().startswith("y"):
            print("Sequences accepted.")
            return True

        elif query.lower().startswith("n"):
            raise ValueError("This system only accepts protein sequences.")

        else:
            raise ValueError("Response unintelligible.")

    else:
            raise ValueError("Only protein sequences accepted.")
