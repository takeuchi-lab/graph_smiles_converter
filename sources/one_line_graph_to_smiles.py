import networkx as nx
import pysmiles
import pandas as pd
# RDKitのインポート
from rdkit import Chem
from rdkit.Chem import rdmolops
from typing import List

# bond記号→order値の変換テーブル
BOND_SYMBOL_TO_ORDER = {'-': 1, '=': 2, ':': 1.5, '#': 3}

# order値→RDKit BondType
ORDER_TO_BONDTYPE = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    1.5: Chem.BondType.AROMATIC,
    3: Chem.BondType.TRIPLE
}

ELEMENT_TO_NUM = Chem.GetPeriodicTable().GetAtomicNumber

def graph_reader(graph_str: str) -> nx.Graph:
    """
    1行形式のグラフ文字列から networkx.Graph を構築する。
    例: "v 0 C v 1 N v 2 C ... e 0 1 - e 1 2 ="
    ノード: "v" の後にノードIDとラベルが続く
    エッジ: "e" の後にソース, ターゲット, bond記号（例 '-', '=', ':'）が続く
    Returns:
        nx.Graph: 構築されたグラフ
    """
    tokens = graph_str.strip().split()
    G = nx.Graph()
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == 'v':
            node_id = int(tokens[i+1])
            label = tokens[i+2]
            G.add_node(node_id, element=label)
            i += 3
        elif token == 'e':
            source = int(tokens[i+1])
            target = int(tokens[i+2])
            bond = tokens[i+3]
            order = BOND_SYMBOL_TO_ORDER.get(bond, 1)
            G.add_edge(source, target, order=order)
            i += 4
        else:
            i += 1
    return G

# networkxグラフ→RDKit Mol
def nx_to_rdkit_mol(G: nx.Graph) -> Chem.Mol:
    """
    networkxグラフをRDKit Molオブジェクトに変換する。
    Args:
        G (nx.Graph): 元のグラフ
    Returns:
        Chem.Mol: RDKit分子オブジェクト
    """
    mol = Chem.RWMol()
    node_to_idx = {}
    # ノード追加
    for n, data in G.nodes(data=True):
        elem = data['element']
        atom = Chem.Atom(elem)
        idx = mol.AddAtom(atom)
        node_to_idx[n] = idx
    # エッジ追加
    for u, v, data in G.edges(data=True):
        order = data.get('order', 1)
        # floatの誤差を吸収
        if abs(order - 1.5) < 0.1:
            bondtype = Chem.BondType.AROMATIC
        elif abs(order - 1) < 0.1:
            bondtype = Chem.BondType.SINGLE
        elif abs(order - 2) < 0.1:
            bondtype = Chem.BondType.DOUBLE
        elif abs(order - 3) < 0.1:
            bondtype = Chem.BondType.TRIPLE
        else:
            bondtype = Chem.BondType.SINGLE
        mol.AddBond(node_to_idx[u], node_to_idx[v], bondtype)
    # 芳香族フラグの自動設定
    rdmolops.Kekulize(mol, clearAromaticFlags=True)
    return mol

def graph_to_smiles(graph: nx.Graph) -> str:
    """
    networkxグラフからSMILES文字列を生成する。
    Args:
        graph (nx.Graph): 入力グラフ
    Returns:
        str: SMILES文字列
    """
    mol = nx_to_rdkit_mol(graph)
    try:
        smiles = Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        smiles = ''
    return smiles

def make_smiles_from_graphs(graphs: pd.Series) -> List[str]:
    """
    複数のグラフ文字列からSMILESリストを生成する。
    Args:
        graphs (pd.Series): グラフ文字列のリスト
    Returns:
        List[str]: SMILES文字列リスト
    """
    smiles_list = []
    for graph in graphs:
        nx_graph = graph_reader(graph)
        smiles_list.append(graph_to_smiles(nx_graph))
    return smiles_list

def make_smiles_from_csv(
    csv_path: str,
    graph_column_name: str,
    other_column_names: List[str],
    output_csv_path: str
) -> None:
    """
    CSVファイルからグラフ列をSMILES列に変換し、他の列とともに新しいCSVとして出力する。
    Args:
        csv_path (str): 入力CSVパス
        graph_column_name (str): グラフ列名
        other_column_names (List[str]): 追加で出力する列名
        output_csv_path (str): 出力CSVパス
    """
    df = pd.read_csv(csv_path)
    smiles_list = make_smiles_from_graphs(df[graph_column_name])
    output_df = pd.DataFrame(columns=["smiles"] + other_column_names)
    output_df["smiles"] = smiles_list
    for column_name in other_column_names:
        output_df[column_name] = df[column_name]
    output_df.to_csv(output_csv_path, index=False)

if __name__ == "__main__":
    # サンプル実行例
    make_smiles_from_csv("Onur_Ebselen_graph.csv", "graph", [], "Onur_Ebselen_smiles.csv")