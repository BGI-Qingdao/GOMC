import sqlite3
import pandas as pd
import numpy as np
from os import path
from sys import argv

def fetch_gcf_models(result_folder):
    print("loading gcf features...")
    with sqlite3.connect(path.join(result_folder, "result/data.db")) as con:
        cur = con.cursor()
        gcf_ids = list(set([row[0] for row in cur.execute("select gcf_id from gcf_membership order by gcf_id asc").fetchall()]))
        hmm_ids = [row[0] for row in cur.execute("select id, name from hmm where db_id=1 order by id asc").fetchall()]
        gcf_features = pd.DataFrame(
            np.zeros((len(gcf_ids), len(hmm_ids)), dtype=np.uint8),
            index=gcf_ids,
            columns=hmm_ids
        )
        for gcf_id, hmm_id, value in cur.execute((
            "select gcf_id, hmm_id, value"
            " from gcf_models,hmm"
            " where hmm.id=gcf_models.hmm_id"
            " and hmm.db_id=1"
        )).fetchall():
            gcf_features.at[gcf_id, hmm_id] = value
    return gcf_features

def get_hmm_names(result_folder):
    with sqlite3.connect(path.join(result_folder, "result/data.db")) as con:
        cur = con.cursor()
        ids, names = list(zip(*cur.execute("select id, name from hmm where db_id=1 order by id asc").fetchall()))
        hmm_names = pd.Series(names, index=ids)
        return hmm_names
    
def main():
    try:
        bigslice_result_folder = argv[1]
        output_csv_path = argv[2]
    except:
        print("usage: python extract_gcf_models_matrix.py <bigslice_result_folder> <output_tsv_path>")
        return 1

    gcf_features = fetch_gcf_models(bigslice_result_folder)
    gcf_features.columns = get_hmm_names(bigslice_result_folder)

    print ("saving to file...")
    gcf_features.to_csv(output_csv_path, sep="\t")
    return 0

if __name__ == "__main__":
    main()