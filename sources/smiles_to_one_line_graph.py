from pysmiles import read_smiles
import pandas as pd
import sys
import csv
import os
import re
from collections import defaultdict
import networkx as nx
import argparse
from typing import Tuple, List

EXTRACTED_NODE_INFO = 'element'
EXTRACTED_EDGE_INFO = 'order'

def dict_order(ls: List[str], init: int = 0) -> dict:
    """
    リストから要素→インデックスの辞書を作成
    """
    return {e: i + init for i, e in enumerate(ls)}

NUMBER2ELEMENT = [None, "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]
ELEMENT2NUMBER = dict_order(NUMBER2ELEMENT)
NUMBER2BOND = {2: '-', 3: ':', 4: '=', 6: '#'}

def iloc_num_name(df: pd.DataFrame, row_num: int, col_name: str):
    """
    DataFrameの行番号・列名から値を取得
    """
    return df.loc[df.index[row_num], col_name]

def convert_graph_to_csv_format(graph: nx.Graph, objective_value) -> Tuple:
    """
    NetworkXグラフをCSV形式の1行文字列に変換
    Returns:
        (目的変数, グラフ文字列)
    """
    vertices = []
    edges = []
    nodes = graph.nodes(data=True)
    for node in nodes:
        node_label = ELEMENT2NUMBER[node[1][EXTRACTED_NODE_INFO]]
        if "hcount" in node[1]:
            node_label += node[1]["hcount"] * 1000
        vertices.append(f"v {node[0]} {NUMBER2ELEMENT[node_label%1000]}")
    edges_data = graph.edges(data=True)
    for edge in edges_data:
        edge_info = edge[2][EXTRACTED_EDGE_INFO]
        if round(edge_info, 2) == 1:
            bond_symbol = '-'
        elif round(edge_info, 2) == 2:
            bond_symbol = '='
        elif round(edge_info, 2) == 1.5:
            bond_symbol = ':'
        elif round(edge_info, 2) == 3:
            bond_symbol = '#'
        else:
            bond_symbol = '-'
        edges.append(f"e {edge[0]} {edge[1]} {bond_symbol}")
    graph_str = ' '.join(vertices + edges)
    return objective_value, graph_str

def smiles_to_csv(
    input_file: str,
    output_file: str,
    hydrogen: int,
    colname_smiles: str,
    objective_column: str
) -> None:
    """
    SMILESファイルをCSV形式に直接変換
    Args:
        input_file: 入力ファイル（CSVまたはExcel）
        output_file: 出力CSVファイル
        hydrogen: 水素を含めるか（1:含める, 0:含めない）
        colname_smiles: SMILES列名
        objective_column: 目的変数列名
    """
    hydrogen = int(hydrogen)
    if input_file[-5:] == '.xlsx':
        data = pd.read_excel(input_file, header=0, index_col=None)
    else:
        data = pd.read_csv(input_file, header=0, index_col=None)
    columns_not_found = []
    if colname_smiles not in data.columns:
        columns_not_found.append(colname_smiles)
    if objective_column not in data.columns:
        columns_not_found.append(objective_column)
    if len(columns_not_found) > 0:
        raise RuntimeError(f'Specified columns ({", ".join(columns_not_found)}) not included in the file "{input_file}"')
    csv_data = []
    for i in range(len(data)):
        smstr = iloc_num_name(data, i, colname_smiles)
        objective_value = iloc_num_name(data, i, objective_column)
        smstr_split = smstr.split(".")
        if len(smstr_split) > 1:
            smstr = max(smstr_split, key=lambda s: len(s))
        graph = read_smiles(smstr, explicit_hydrogen=hydrogen)
        obj_val, graph_str = convert_graph_to_csv_format(graph, objective_value)
        csv_data.append((obj_val, graph_str))
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(["objective", "graph"])
        for objective_val, graph_str in csv_data:
            writer.writerow([objective_val, graph_str])

def main() -> None:
    """
    コマンドライン引数をパースし、SMILES→1行グラフCSV変換を実行する。
    """
    parser = argparse.ArgumentParser(description="Convert SMILES file directly to CSV format containing multiple molecular structures.")
    parser.add_argument("input_file", help="Path to the input SMILES file (CSV or Excel)")
    parser.add_argument("hydrogen", help="Include hydrogen (1) or not (0)")
    parser.add_argument("smiles_column", help="Name of the column containing SMILES")
    parser.add_argument("objective_column", help="Name of the objective column")
    parser.add_argument("output_file", help="Path to the output CSV file")
    args = parser.parse_args()
    smiles_to_csv(args.input_file, args.output_file, args.hydrogen, args.smiles_column, args.objective_column)
    print(f"Conversion completed: {args.input_file} -> {args.output_file}")

if __name__ == '__main__':
    if len(sys.argv) < 6:
        sys.stderr.write("Usage: smiles_to_csv [SMILES_TABLE_FILE] [INCLUDE_HYDROGEN] [SMILES_COLUMN] [OBJECTIVE_COLUMN] [OUTPUT_FILE]\n")
        sys.exit(-1)
    main() 