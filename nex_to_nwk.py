import dendropy
import argparse

def transfer_posterior(tree):
    for node in tree.preorder_node_iter():
        if not node.is_leaf():
            posterior = node.annotations.get_value("posterior")
            if posterior is not None:
                node.label = str(posterior)

def main():
    parser = argparse.ArgumentParser(
        description="Convert a BEAST-generated Nexus tree to Newick format with an option to include posterior probabilities as support values."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input Nexus tree file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output Newick tree file")
    parser.add_argument("--include-posterior", action="store_true",
                        help="If set, transfer posterior probabilities from annotations to internal node labels")
    args = parser.parse_args()

    try:
        tree = dendropy.Tree.get(path=args.input, schema="nexus")
    except Exception as e:
        print("Error reading Nexus file:", e)
        return

    if args.include_posterior:
        transfer_posterior(tree)

    newick_str = tree.as_string(schema="newick")

    with open(args.output, "w") as out_file:
        out_file.write(newick_str)

    print(f"Conversion successful! Newick tree saved as: {args.output}")

if __name__ == '__main__':
    main()
