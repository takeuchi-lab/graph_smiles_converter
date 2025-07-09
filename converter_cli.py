"""
任意の化合物ファイル（SMILES, 1行グラフ, gSpan形式）を相互変換する統合CLI

- SMILES→1行グラフCSV
- 1行グラフCSV→gSpan形式
- gSpan形式→1行グラフCSV
- 1行グラフCSV→SMILES

Usage例:
  python converter_cli.py --mode smiles_to_graph --input input.csv --output output.csv --smiles-column smiles --objective-column y --hydrogen 1
  python converter_cli.py --mode graph_to_gspan --input input.csv --output output.graph
  python converter_cli.py --mode gspan_to_graph --input input.graph --output output.csv
  python converter_cli.py --mode graph_to_smiles --input input.csv --output output.csv
"""
import argparse
import sys
from typing import List

from sources.smiles_to_one_line_graph import smiles_to_csv
from sources.one_line_graph_to_gspan_format import convert_csv_to_graph
from sources.gspan_format_to_one_line_graph import convert_gspan_to_csv
from sources.one_line_graph_to_smiles import make_smiles_from_csv

def main():
    parser = argparse.ArgumentParser(
        description="SMILES, 1行グラフ, gSpan形式の相互変換ツール"
    )
    parser.add_argument('--mode', required=True, choices=['smiles_to_graph', 'graph_to_gspan', 'gspan_to_graph', 'graph_to_smiles'],
                        help='変換モード: smiles_to_graph, graph_to_gspan, gspan_to_graph, graph_to_smiles')
    parser.add_argument('--input', required=True, help='入力ファイルパス')
    parser.add_argument('--output', required=True, help='出力ファイルパス')
    # SMILES→1行グラフ用
    parser.add_argument('--smiles-column', default='smiles', help='SMILES列名 (smiles_to_graph用)')
    parser.add_argument('--objective-column', default='objective', help='目的変数列名 (smiles_to_graph用)')
    parser.add_argument('--hydrogen', default='0', help='水素を含めるか (1:含める, 0:含めない) (smiles_to_graph用)')
    # 1行グラフ→gSpan用
    parser.add_argument('--graph-column', default='graph', help='グラフ列名 (graph_to_gspan, graph_to_smiles用)')
    parser.add_argument('--spacing', action='store_true', help='gSpan出力時にグラフ間に空行を挿入')
    parser.add_argument('--target-column', default=None, help='gSpan出力時にヘッダーへ目的変数（ターゲット値）を出力する列名（例: objective, y など, オプション）')
    # gSpan→1行グラフ用
    parser.add_argument('--gspan-target-column', default='objective', help='gSpan入力時の目的変数列名 (gspan_to_graph用)')
    # 1行グラフ→SMILES用
    parser.add_argument('--other-columns', nargs='*', default=[], help='SMILES出力時に一緒に出力する他の列名 (graph_to_smiles用)')

    args = parser.parse_args()

    if args.mode == 'smiles_to_graph':
        smiles_to_csv(
            input_file=args.input,
            output_file=args.output,
            hydrogen=args.hydrogen,
            colname_smiles=args.smiles_column,
            objective_column=args.objective_column
        )
        print(f"[OK] SMILES→1行グラフCSV変換完了: {args.input} → {args.output}")

    elif args.mode == 'graph_to_gspan':
        convert_csv_to_graph(
            input_file=args.input,
            output_file=args.output,
            graph_column=args.graph_column,
            add_spacing=args.spacing,
            target_column=args.target_column
        )
        print(f"[OK] 1行グラフCSV→gSpan変換完了: {args.input} → {args.output}")

    elif args.mode == 'gspan_to_graph':
        convert_gspan_to_csv(
            input_file=args.input,
            output_file=args.output,
            target_column_name=args.gspan_target_column
        )
        print(f"[OK] gSpan→1行グラフCSV変換完了: {args.input} → {args.output}")

    elif args.mode == 'graph_to_smiles':
        make_smiles_from_csv(
            csv_path=args.input,
            graph_column_name=args.graph_column,
            other_column_names=args.other_columns,
            output_csv_path=args.output
        )
        print(f"[OK] 1行グラフCSV→SMILES変換完了: {args.input} → {args.output}")

    else:
        print('未対応のmodeです', file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main() 