import ui_script as uis
from Bio import SeqIO
import programme_settings as ps
import blast_101_search as blast
from ui_script import argument_taker


def analysis(algo, database_path, sequence_path):
    if algo == "blast101":
        print("Running BLAST101..")
        if "DEFAULT" in ps.settings and ps.settings["DEFAULT"].get("database") is not None:
            ps.settings["DEFAULT"]["database"] = database_path

        if ("DEFAULT" in ps.settings and "query_sequence" in ps.settings["DEFAULT"]) is not None:
            ps.settings["DEFAULT"]["query_sequence"] = sequence_path

        blast.blast101_run()
        print(f"Executed BLAST101.")

    elif algo == "smith-waterman":
        print("Running Smith-Waterman...")
        if "DEFAULT" in ps.settings and ps.settings["DEFAULT"].get("database") is not None:
            ps.settings["DEFAULT"]["database"] = database_path

        if ("DEFAULT" in ps.settings and "query_sequence" in ps.settings["DEFAULT"]) is not None:
            ps.settings["DEFAULT"]["query_sequence"] = sequence_path
        import smith_waterman_search as sws
        sws.processSW()
        print("Executed Smith-Waterman...")

    elif algo == "statistical_analysis":
        print("Conducting statistical analysis...")
        if "DEFAULT" in ps.settings and ps.settings["DEFAULT"].get("database") is not None:
            ps.settings["DEFAULT"]["database"] = database_path

        if ("DEFAULT" in ps.settings and "query_sequence" in ps.settings["DEFAULT"]) is not None:
            ps.settings["DEFAULT"]["query_sequence"] = sequence_path

        import calc_bit_and_evalues as cbe
        cbe.build_fit()
        cbe.test()
        print("Finished statistical analysis.")

    return

def ui_structure():
    #Step 1: Take input from the user
    args = argument_taker()
    database_path = args.database
    sequence_path = args.sequence
    algo = args.method
    run_tests2 = args.tests

    #Step 2: Run the method
    analysis(algo, database_path, sequence_path)
    #Step 3 (Optional): Allow the user to test the commandline
    if run_tests2 == "true":
        from test_script import waitt
        waitt()
    return



if __name__ == "__main__":
    ui_structure()
