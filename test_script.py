################################################
# Part 2: Testing the C #
################################################

import os, builtins
from ui_script import does_it_exist, but_does_it_fasta, protein_or_not
from pleasework import analysis
import blast_101_search as blast
import smith_waterman_p as swp

#Test files
functional_db = "tempdb.fasta"
functional_query = "tempq.fasta"
with open(functional_db, 'w') as fdb:
    fdb.write(f">seq1\nMALWMRLLPLLALLALWGPDPAAAEQKNAKDEQMQKALEDGQARFDLLESLVD\n>seq2\nMVKVGVNGFGRIGRLVTRAAFNSGKVDIVLDSGDGVTHVVAGKHAYNVHAYLNS\n>seq3\nMFKTFVLLLQLALRSLPAGKLRLLTAAKPDHQGLLRSLGVAGIVAGKMPVNGE")

with open(functional_query, 'w') as fq:
    fq.write(f">query1\nMVKVGVNGFGRIGRLVTRAAFNSGKVDIVLDSGDGVTHVVAGKHAYNVHAYLNS\n>query2\nMFKTFVLLLQLALRSLPAGKLRLLTAAKPDHQGLLRSLGVAGIVAGKMPVNGE")

empty_fasta = "empty.fasta"
open(empty_fasta, 'w').close()

not_fasta = "notafasta.txt"
with open(not_fasta, 'w') as nf:
    nf.write(f">query1\nMVKVGVNGFGRIGRLVTRAAFNSGKVDIVLDSGDGVTHVVAGKHAYNVHAYLNS\n>query2\nMFKTFVLLLQLALRSLPAGKLRLLTAAKPDHQGLLRSLGVAGIVAGKMPVNGE")

not_fasta2 = "alsonotafasta.txt"
with open(not_fasta2, 'w') as nf2:
    nf2.write(f">seq1\nMALWMRLLPLLALLALWGPDPAAAEQKNAKDEQMQKALEDGQARFDLLESLVD\n>seq2\nMVKVGVNGFGRIGRLVTRAAFNSGKVDIVLDSGDGVTHVVAGKHAYNVHAYLNS\n>seq3\nMFKTFVLLLQLALRSLPAGKLRLLTAAKPDHQGLLRSLGVAGIVAGKMPVNGE")

no_file = " "

fake_dna = "fakedna.fasta"
with open(fake_dna, 'w') as fd:
    fd.write(f">seq1\nATGCATGCATGCATGCATGC\nseq2\nTGCAATGCATGCATGCATG\n>seq3\nATGTGCAAGCTAGCATGCAT")

dna = "dna.fasta"
with open(dna, 'w') as fd2:
    fd2.write(">query1\nATGCATGCATGC\n>query2\nTGCAATGCATGC\n")

weird_fasta = "weird.fasta"
with open(weird_fasta, 'w') as wf:
    wf.write(">seq1\nMKTAYIAKQR#ISFVKSHFSRQDQQGAMVIGRVPAGAVVGMVGQETQDAAGGAYKVC\n>seq2\nMAYIVTGGVLSAGF*LLVLVGGSGGAAGVVSPITPAGTVVSHKAKANNSFPFRLFGT\n>seq3\nMTGGVILASGTGKAKNRMFGSRSIT^GNVFSPFGFVAVSGAGDHRFQKAAEKLNVNS")






###############################
# Section 1: Input validation #
###############################

#This tests that the does_it_exist() function works properly
def exist_inator():
    #The scenario where both files exist and are populated
    try:
        does_it_exist(functional_db, functional_query)
        print("PASS: Files exist and are populated.")
    except Exception as e:
        print(f"FAIL: Files rejected -> {str(e)}")

    #The scenario where at least one of the files doesn't exist
    try:
        does_it_exist(functional_db, no_file)
    except FileNotFoundError as fnf:
        print(f"PASS: {str(fnf)}")

    #The scenario where both files exist but one is empty
    try:
        does_it_exist(functional_db, empty_fasta)
    except ValueError as ve:
        print(f"PASS: {str(ve)}")

    return

#This tests if the but_does_it_fasta() function works properly
def fasta_inator():
    #The scenario where both files are FASTA files
    try:
        but_does_it_fasta(functional_db, functional_query)
        print("PASS: Files exist and are populated.")
    except Exception as e:
        print(f"FAIL: Files rejected -> {str(e)}")

    # the scenario where one of the files id not a FASTA files
    try:
        but_does_it_fasta(functional_db, not_fasta)
    except ValueError as ve:
        print(f"PASS: {str(ve)}")

    #the scenario where both of the files are not FASTA files
    try:
        but_does_it_fasta(not_fasta2, not_fasta)
    except ValueError as ve:
        print(f"PASS: {str(ve)}")

    return
original_input = builtins.input

#This tests if the protein_or_not() function works properly
def proto_inator():
    #The scenario where both sequences are protein sequences
    try:
        protein_or_not(functional_db, functional_query)#
        print("PASS: Protein sequences accepted.")
    except Exception:
        print("FAIL: Sequences rejected.")

    #The scenario where one of the sequences is a dna sequence
    #When the programme asks us if we're sure about our sequences, we say no.
    builtins.input = lambda _: "n"
    try:
        protein_or_not(functional_db, dna)
    except ValueError as ve:
        print(f"PASS: {str(ve)}")

    #The scenario where one of the sequences contains a weird character
    try:
        protein_or_not(functional_db, weird_fasta)
    except ValueError:
        print("PASS: One of these files contains a non-amino-acid character.")

    return


################################
# Section 2: System Validation #
################################
#This function ensures that the pipeline works properly
def schmipeline_inator():
    #Create temporary files to use as input

    #When the CLI asks us if we wish to proceed with our files, we say yes.
    builtins.input = lambda _: "y"

    try:
        #Testing out BLAST
        print("Testing the BLAST method...")
        analysis("blast101", functional_db, functional_query)
        print("BLAST101 Executed Successfully")
    except Exception as e:
        print(f"FAIL: Pipeline crashed -> {type(e).__name__}: {str(e)}")

    builtins.input = original_input
    question = input("Are you sure you would like to test Smith-Waterman (y/n)? This may take a very long time.")
    if question.lower().startswith("y"):
        try:
        #Testing out Smith-Waterman
            print("Testing the Smith-Waterman method...")
            analysis("smith-waterman", functional_db, functional_query)
            print("Smith Waterman Executed Successfully")
        except Exception as e:
            print(f"FAIL: Pipeline crashed -> {type(e).__name__}: {str(e)}")
    else:
        print("Skipping Smith-Waterman...")

    try:
        #Testing out statistical analysis
        print("Testing the statistical analyses...")
        analysis("statistical_analysis", functional_db, functional_query)
        print("Statistical Analyses Executed Successfully")

    except Exception as e:
        print(f"FAIL: Pipeline crashed -> {type(e).__name__}: {str(e)}")

    return

##############################
# Section 2B: BLAST Specific #
##############################
#Testing the score verification in the extend_diagonal() function from blast_101_search
def scorolo():
    starting_pos = (0, 0)

    #Test 1: Perfect Match
    try:
        blast.word_size = 2
        s0 = "AC"
        s1 = "AC"
        score1 = blast.extend_diagonal(starting_pos, s0, s1)
        if score1 == 13:
            print(f"Perfect Alignment Test PASSED. BLAST Alignment Score: {score1}.")
        else:
            print(f"Perfect Alignment Test FAILED. BLAST Alignment Score: {score1} ")
    except Exception as e:
        print(f"Crash: Incorrect BLAST Alignment Score -> {type(e).__name__}: {str(e)}")

    #Test 2: Complete Mismatch
    try:
        s0 = "AC"
        s1 = "DE"
        score2 = blast.extend_diagonal(starting_pos, s0, s1)
        if score2 <= 0:
            print(f"Mismatched Sequences Test PASSED. BLAST Score: {score2}. Expected Score: <= 0")
        else:
            print(f"Mismatched Sequences Test FAILED. BLAST Score: {score2}. Expected Score: <= 0")
    except Exception:
        print(f"FAIL: Mismatch test crashed")

    #Test 3: Weak Similarity
    try:
        s0 = "AC"
        s1 = "AD"
        score3 = blast.extend_diagonal(starting_pos, s0, s1)
        if 0 < score2 < 5:
            print(f"Weak Similarity Test PASSED. BLAST Score: {score2}. Excpected Score: 0 < score < 5")
        else:
            print(f"Weak Similarity Test FAILED. BLAST Score: {score2}. Excpected Score: 0 < score < 5")
    except Exception as e:
        print(f"CRASH: Weak Similarity Test -> {type(e).__name__}: {str(e)}")

    #Test 4: The Biologically Implausible
    try:
        s0 = "KR"
        s1 = "VL"
        score4 = blast.extend_diagonal(starting_pos, s0, s1)
        if score2 < 0:
            print(f"Biological Implausibility Test PASSED. BLAST Score: {score2}. Expected Score: < 0")
        else:
            print(f"Biological Implausibility Test FAILED. BLAST Score: {score2}. Expected Score: < 0")
    except Exception as e:
        print(f"CRASH: Implausible Pair Test -> {type(e).__name__}: {str(e)}")

    #Test 5: Mismatched Lengths
    try:
        s0 = "ACD"
        s1 = "AC"
        score5 = blast.extend_diagonal(starting_pos, s0, s1)
        if score5 == 13:
            print(f"Mismatched Length Test PASSED. BLAST Score: {score5}. Expected Score: 13")
        else:
            print(f"Mismatched Length Test FAILED. BLAST Score: {score5}. Expected Score: 13")
    except Exception as e:
        print(f"CRASH: Mismatch Length Test -> {type(e).__name__}: {str(e)}")

    #Mismatched Word Length
    # Test 5: Mismatched Lengths
    try:
        blast.word_size = 1
        blast.tin_test = 14
        s0 = "ACC"
        s1 = "ACC"
        score6 = blast.extend_diagonal(starting_pos, s0, s1)
        if score6 < blast.tin_test:
            print(f"Mismatched Word Length Test PASSED. BLAST Score: {score6}. Expected Score: 13")
        else:
            print(f"Mismatched Word Length Test FAILED. BLAST Score: {score6}. Expected Score: 13")
    except Exception as e:
        print(f"CRASH: Mismatch Word Length Test -> {type(e).__name__}: {str(e)}")
    return

#######################################
# Section 2B: Smith-Waterman Specific #
#######################################
def sw_tester():
    #Test 1: Perfect Match
    try:
        s1 = "AC"
        s2 = "AC"
        score1 = swp.perform_smith_waterman(s1, s2, print_a = True)
        if score1 > 0:
            print(f"Perfect Alignment Test PASSED. SW Alignment Score: {score1}.")
        else:
            print(f"Perfect Alignment Test FAILED. SW Alignment Score: {score1}.")
    except Exception as e:
        print(f"Crash: Incorrect SW Alignment Score -> {type(e).__name__}: {str(e)}")

    #Test 2: Complete Mismatch
    try:
        s1 = "AC"
        s2 = "DE"
        score2 = swp.perform_smith_waterman(s1, s2, print_a = True)
        if score2 < 0:
            print(f"Complete Mismatch Test PASSED. SW Alignment Score: {score2}.")
        else:
            print(f"Complete Mismatch Test FAILED. SW Alignment Score: {score2}.")
    except Exception as e:
        print(f"Crash: Incorrect SW Alignment Score -> {type(e).__name__}: {str(e)}")

    #Test 3: Weak Similarity
    try:
        s1 = "AC"
        s2 = "AD"
        score3 = swp.perform_smith_waterman(s1, s2, print_a = True)
        if 0 < score3 < 5:
            print(f"Weak Similarity Test PASSED. SW Alignment Score: {score3}.")
        else:
            print(f"Weak Similarity Test FAILED. SW Alignment Score: {score3}.")
    except Exception as e:
        print(f"Crash: Incorrect SW Alignment Score -> {type(e).__name__}: {str(e)}")

    # Test 4: Gap Penalty
    try:
        s1 = "A-C"
        s2 = "ACC"
        score4 = swp.perform_smith_waterman(s1, s2, print_a=True)
        if score4 == 9:
            print(f"Gap Penalty Test PASSED. SW Alignment Score: {score4}.")
        else:
            print(f"Gap Penalty Test FAILED. SW Alignment Score: {score4}.")
    except Exception as e:
        print(f"Crash: Incorrect SW Alignment Score -> {type(e).__name__}: {str(e)}")
    return

#########################
# The Tests, Altogether #
#########################
def waitt():
    print("Running tests...")
    print("Phase 1: Testing the Command Line Interface...")
    exist_inator()
    fasta_inator()
    proto_inator()
    print("Phase 2: Testing the Pipeline")
    schmipeline_inator()
    print("Phase 3: Testing BLAST Functions...")
    scorolo()
    print("Phase 4: Testing Smith-Waterman Functions...")
    sw_tester()

    # Delete the temporary files
    os.remove(functional_query)
    os.remove(functional_db)
    os.remove(empty_fasta)
    os.remove(fake_dna)
    os.remove(dna)
    os.remove(weird_fasta)
    os.remove(not_fasta)
    os.remove(not_fasta2)
    return