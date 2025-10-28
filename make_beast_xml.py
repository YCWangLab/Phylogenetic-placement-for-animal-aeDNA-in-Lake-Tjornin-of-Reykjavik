import argparse
from Bio import SeqIO
from xml.etree import ElementTree as ET
import xml.dom.minidom as minidom
import re


def load_xml_template(xml_file):
    with open(xml_file, "r", encoding="utf-8") as f:
        xml_content = f.read()

    match = re.match(r'<\?xml.*?\?>', xml_content)
    xml_declaration = match.group(0) if match else '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'

    tree = ET.ElementTree(ET.fromstring(xml_content))
    root = tree.getroot()

    return tree, root, xml_declaration


def replace_prefix_in_xml(root, old_prefix, new_prefix):
    for elem in root.iter():
        for key, value in elem.attrib.items():
            if old_prefix in value:
                new_value = re.sub(rf'([:@.]?){old_prefix}', rf'\1{new_prefix}', value)
                elem.set(key, new_value)


def update_mcmc_params(root, chain_length, store_every):
    run_elem = root.find(".//run[@id='mcmc']")
    if run_elem is not None:
        run_elem.set("chainLength", str(chain_length))
        run_elem.set("storeEvery", str(store_every))
        print(f"Updated MCMC params: chainLength={chain_length}, storeEvery={store_every}")
    else:
        print("Warning: <run id='mcmc'> element not found.")


def insert_fasta_into_xml(fasta_file, root, alignment_id):
    data_block = root.find(f".//data[@id='{alignment_id}']")
    if data_block is None:
        print(f"Warning: <data id='{alignment_id}'> not found.")
        return

    # Clear existing sequences
    for child in list(data_block):
        data_block.remove(child)

    # Insert new sequences from the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    for seq in sequences:
        sequence_elem = ET.SubElement(data_block, "sequence")
        sequence_elem.set("id", f"seq_{seq.id}")
        sequence_elem.set("spec", "Sequence")
        sequence_elem.set("taxon", seq.id)
        sequence_elem.set("totalcount", "4")
        sequence_elem.set("value", str(seq.seq))

    print(f"Inserted {len(sequences)} sequences into <data id='{alignment_id}'>.")


def save_xml(tree, output_file, xml_declaration):
    xml_str = ET.tostring(tree.getroot(), encoding="UTF-8")

    parsed = minidom.parseString(xml_str)
    formatted_xml = parsed.toprettyxml(indent="  ")

    if xml_declaration:
        formatted_xml = formatted_xml.split('\n', 1)[1]  # Remove auto-generated declaration
        formatted_xml = f"{xml_declaration}\n{formatted_xml}"

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(formatted_xml)

    print(f"XML file saved as: {output_file}")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Modify BEAST XML files by inserting FASTA sequences and updating MCMC parameters.")
    
    # Define command-line arguments
    parser.add_argument("--xml_template", required=True, help="Path to the input XML template file.")
    parser.add_argument("--fasta_file", required=True, help="Path to the aligned FASTA file.")
    parser.add_argument("--output_file", required=True, help="Path to save the updated XML file.")
    parser.add_argument("--prefix", required=True, help="Prefix to replace 'mito211' in the XML.")
    parser.add_argument("--chain_length", type=int, default=20000000, help="MCMC chain length (default: 20000000).")
    parser.add_argument("--store_every", type=int, default=5000, help="Store interval for MCMC logging (default: 5000).")

    args = parser.parse_args()

    # Execute XML processing steps
    tree, root, xml_declaration = load_xml_template(args.xml_template)

    replace_prefix_in_xml(root, "ultrametric", args.prefix)

    update_mcmc_params(root, args.chain_length, args.store_every)

    insert_fasta_into_xml(args.fasta_file, root, args.prefix)

    save_xml(tree, args.output_file, xml_declaration)


if __name__ == "__main__":
    main()