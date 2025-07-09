import csv
import argparse
import io
from typing import List, Tuple, Optional

# 原子番号から元素記号への変換テーブル
NUMBER2ELEMENT = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
    21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
    41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd",
    61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm",
    101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds",
    111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
}

# 結合数値から記号への変換テーブル
NUMBER2BOND = {2: '-', 3: ':', 4: '=', 6: '#'}

def convert_number_to_element(atomic_number: int) -> str:
    """
    原子番号を元素記号に変換する。
    """
    if atomic_number in NUMBER2ELEMENT:
        return NUMBER2ELEMENT[atomic_number]
    else:
        raise ValueError(f"Unknown atomic number: {atomic_number}")

def convert_number_to_bond(bond_number: int) -> str:
    """
    結合数値を記号に変換する。
    """
    if bond_number in NUMBER2BOND:
        return NUMBER2BOND[bond_number]
    else:
        # デフォルトは単結合
        return '-'

def parse_gspan_graph(lines: List[str]) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str, str]], Optional[str], Optional[str]]:
    """
    gSpan形式のグラフ行を解析し、頂点と辺の情報を抽出する。
    Returns:
        vertices: List of (vertex_id, element)
        edges: List of (from_id, to_id, bond)
        target_name: 目的変数の列名（存在する場合）
        target_value: 目的変数の値（存在する場合）
    """
    vertices = []
    edges = []
    target_name = None
    target_value = None
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        parts = line.split()
        if parts[0] == 't':
            # ヘッダー行: t # 0 [target_name] [target_value]
            if len(parts) >= 5:
                target_name = parts[3]
                target_value = parts[4]
        elif parts[0] == 'v':
            # 頂点行: v vertex_id atomic_number
            vertex_id = parts[1]
            atomic_number = int(parts[2])
            element = convert_number_to_element(atomic_number)
            vertices.append((vertex_id, element))
        elif parts[0] == 'e':
            # 辺行: e from_id to_id bond_number
            from_id = parts[1]
            to_id = parts[2]
            bond_number = int(parts[3])
            bond = convert_number_to_bond(bond_number)
            edges.append((from_id, to_id, bond))
    
    return vertices, edges, target_name, target_value

def graph_to_csv_format(vertices: List[Tuple[str, str]], edges: List[Tuple[str, str, str]]) -> str:
    """
    頂点と辺の情報を1行グラフCSV形式の文字列に変換する。
    Returns:
        1行グラフ形式の文字列
    """
    vertex_strs = [f"v {vid} {element}" for vid, element in vertices]
    edge_strs = [f"e {from_id} {to_id} {bond}" for from_id, to_id, bond in edges]
    return ' '.join(vertex_strs + edge_strs)

def convert_gspan_to_csv(
    input_file: str,
    output_file: str,
    target_column_name: str = "objective"
) -> None:
    """
    gSpan形式のグラフファイルを1行グラフCSVファイルに変換する。
    Args:
        input_file: 入力gSpanファイルパス
        output_file: 出力CSVファイルパス
        target_column_name: 目的変数の列名（デフォルト: objective）
    """
    graphs_data = []
    current_graph_lines = []
    
    # gSpanファイルを読み込み、グラフごとに分割
    with open(input_file, 'r', encoding='utf-8') as graphfile:
        for line in graphfile:
            line = line.strip()
            if line.startswith('t #') and current_graph_lines:
                # 新しいグラフの開始
                graphs_data.append(current_graph_lines)
                current_graph_lines = []
            current_graph_lines.append(line)
        
        # 最後のグラフを追加
        if current_graph_lines:
            graphs_data.append(current_graph_lines)
    
    # CSVファイルに出力（LF改行コードで固定）
    with open(output_file, 'w', newline='\n', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, lineterminator='\n')
        
        # ヘッダー行を決定
        has_target = any(
            len(lines) > 0 and len(lines[0].split()) >= 5 
            for lines in graphs_data
        )
        
        if has_target:
            writer.writerow([target_column_name, "graph"])
        else:
            writer.writerow(["graph"])
        
        # 各グラフを処理
        for graph_lines in graphs_data:
            vertices, edges, target_name, target_value = parse_gspan_graph(graph_lines)
            graph_str = graph_to_csv_format(vertices, edges)
            
            if has_target and target_value is not None:
                writer.writerow([target_value, graph_str])
            else:
                writer.writerow([graph_str])

def main() -> None:
    """
    コマンドライン引数をパースし、gSpan→CSV変換を実行する。
    """
    parser = argparse.ArgumentParser(
        description="Convert .graph (gSpan) file to .csv format containing multiple molecular structures."
    )
    parser.add_argument("input_file", help="Path to the input .graph file")
    parser.add_argument("output_file", help="Path to the output .csv file")
    parser.add_argument("--target-column", default="objective", help="Name of the target/objective column (default: objective)")
    args = parser.parse_args()
    
    convert_gspan_to_csv(args.input_file, args.output_file, args.target_column)
    print(f"Conversion completed: {args.input_file} -> {args.output_file}")

if __name__ == '__main__':
    main() 