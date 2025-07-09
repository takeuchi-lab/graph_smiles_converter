import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import argparse

# コマンドライン引数: csvファイル, 出力ディレクトリ, smiles列名（オプション）
parser = argparse.ArgumentParser(description='SMILESからSVGを生成')
parser.add_argument('csv_path', help='入力CSVファイル')
parser.add_argument('out_dir', help='出力ディレクトリ')
parser.add_argument('--smiles_col', default='smiles', help='SMILESが格納された列名（デフォルト: smiles）')
args = parser.parse_args()

csv_path = args.csv_path
out_dir = args.out_dir
smiles_col = args.smiles_col

# ディレクトリがなければ作成
os.makedirs(out_dir, exist_ok=True)

# CSV読み込み
try:
    df = pd.read_csv(csv_path)
except Exception as e:
    print(f"CSVファイルの読み込みに失敗しました: {e}")
    sys.exit(1)

if smiles_col not in df.columns:
    print(f"CSVに'{smiles_col}'列がありません")
    sys.exit(1)

for idx, row in df.iterrows():
    smiles = row[smiles_col]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"{idx}行目: SMILESのパースに失敗: {smiles}")
        continue
    # SVG描画
    try:
        svg = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(300,300), useSVG=True)
        svg_str = svg if isinstance(svg, str) else svg.data
        out_path = os.path.join(out_dir, f"{idx}.svg")
        with open(out_path, 'w') as f:
            f.write(svg_str)
    except Exception as e:
        print(f"{idx}行目: SVG保存に失敗: {e}") 