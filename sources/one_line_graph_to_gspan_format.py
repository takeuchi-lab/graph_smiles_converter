import csv
import argparse
from typing import List, Tuple

# 元素記号から原子番号への変換テーブル
ELEMENT2NUMBER = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
    "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
    "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
}

# 結合記号から数値への変換テーブル
default_bond_number = 1
BOND2NUMBER = {'-': 2, ':': 3, '=': 4, '#': 6}

def parse_csv_graph_string(graph_str: str) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str, str]]]:
    """
    CSVの1行形式のグラフ文字列を解析し、頂点と辺の情報を抽出する。
    例: 'v 0 C v 1 N e 0 1 -'
    Returns:
        vertices: List of (vertex_id, element)
        edges: List of (from_id, to_id, bond)
    """
    vertices = []
    edges = []
    parts = graph_str.split()
    i = 0
    while i < len(parts):
        if parts[i] == 'v':
            vertex_id = parts[i + 1]
            element = parts[i + 2]
            vertices.append((vertex_id, element))
            i += 3
        elif parts[i] == 'e':
            from_id = parts[i + 1]
            to_id = parts[i + 2]
            bond = parts[i + 3]
            edges.append((from_id, to_id, bond))
            i += 4
        else:
            i += 1
    return vertices, edges

def convert_element_to_number(element: str) -> int:
    """
    元素記号を原子番号に変換する。
    """
    if element in ELEMENT2NUMBER:
        return ELEMENT2NUMBER[element]
    try:
        return int(element)
    except ValueError:
        raise ValueError(f"Unknown element: {element}")

def convert_bond_to_number(bond: str) -> int:
    """
    結合記号を数値に変換する。
    """
    if bond in BOND2NUMBER:
        return BOND2NUMBER[bond]
    try:
        return int(bond)
    except ValueError:
        raise ValueError(f"Unknown bond type: {bond}")

def csv_to_graph_format(vertices: List[Tuple[str, str]], edges: List[Tuple[str, str, str]], graph_id: int, target_name: str = None, target_value: str = None) -> List[str]:
    """
    頂点と辺の情報をgSpan形式のグラフ表現に変換する。
    target_name, target_valueが指定された場合はヘッダーに付加する。
    Returns:
        gSpan形式の各行（リスト）
    """
    if target_name is not None and target_value is not None:
        lines = [f"t # {graph_id} {target_name} {target_value}"]
    else:
        lines = [f"t # {graph_id}"]
    for vertex_id, element in vertices:
        atomic_number = convert_element_to_number(element)
        lines.append(f"v {vertex_id} {atomic_number}")
    for from_id, to_id, bond in edges:
        bond_number = convert_bond_to_number(bond)
        lines.append(f"e {from_id} {to_id} {bond_number}")
    return lines

def convert_csv_to_graph(
    input_file: str,
    output_file: str,
    graph_column: str = "graph",
    add_spacing: bool = False,
    target_column: str = None
) -> None:
    """
    CSVファイルをgSpan形式のグラフファイルに変換する。
    Args:
        input_file: 入力CSVファイルパス
        output_file: 出力グラフファイルパス
        graph_column: グラフデータの列名
        add_spacing: グラフ間に空行を挿入するか
        target_column: 目的変数の列名（指定時はヘッダーに付加）
    """
    with open(input_file, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        if graph_column not in reader.fieldnames:
            available_columns = ", ".join(reader.fieldnames)
            raise ValueError(f"Column '{graph_column}' not found in CSV. Available columns: {available_columns}")
        if target_column is not None and target_column not in reader.fieldnames:
            available_columns = ", ".join(reader.fieldnames)
            raise ValueError(f"Target column '{target_column}' not found in CSV. Available columns: {available_columns}")
        rows = list(reader)
        with open(output_file, 'w', encoding='utf-8') as graphfile:
            for i, row in enumerate(rows):
                graph_str = row[graph_column]
                target_value = row[target_column] if target_column is not None else None
                vertices, edges = parse_csv_graph_string(graph_str)
                graph_lines = csv_to_graph_format(vertices, edges, i, target_column, target_value)
                for line in graph_lines:
                    graphfile.write(line + '\n')
                if add_spacing and i < len(rows) - 1:
                    graphfile.write('\n')

def main() -> None:
    """
    コマンドライン引数をパースし、CSV→gSpan変換を実行する。
    """
    parser = argparse.ArgumentParser(
        description="Convert .csv file back to .graph (gSpan) format containing multiple molecular structures."
    )
    parser.add_argument("input_file", help="Path to the input .csv file")
    parser.add_argument("output_file", help="Path to the output .graph file")
    parser.add_argument("--graph-column", default="graph", help="Name of the column containing graph data (default: graph)")
    parser.add_argument("--spacing", action="store_true", help="Add blank line between each graph data")
    parser.add_argument("--target-column", default=None, help="Name of the column containing target/objective value (optional)")
    args = parser.parse_args()
    convert_csv_to_graph(args.input_file, args.output_file, args.graph_column, args.spacing, args.target_column)
    print(f"Conversion completed: {args.input_file} -> {args.output_file}")

if __name__ == '__main__':
    main() 