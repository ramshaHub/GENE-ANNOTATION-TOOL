# GENE-ANNOTATION-TOOL
import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

def read_fasta(file_path):
    """Read a genomic sequence from a FASTA file."""
    with open(file_path, "r") as file:
        return next(SeqIO.parse(file, "fasta"))

def find_orfs(sequence, min_length=100):
    """Identify open reading frames (ORFs) in the sequence."""
    orfs = []
    start_codon = re.compile('ATG')
    stop_codons = re.compile('TAA|TAG|TGA')
    
    for match in start_codon.finditer(str(sequence)):
        start = match.start()
        for stop in stop_codons.finditer(str(sequence), start):
            if (stop.start() - start) % 3 == 0:
                orf = sequence[start:stop.start()+3]
                if len(orf) >= min_length:
                    orfs.append((start, stop.start()+3, orf))
                break
    return orfs

def annotate_gene(gene_start, gene_end, orf_index):
    """Annotate a gene with its subfeatures."""
    gene_id = f'gene{orf_index+1}'
    subfeatures = [
        {
            'Type': 'exon',
            'Start': gene_start,
            'End': gene_end,
            'ID': f'exon{orf_index+1}',
            'Attributes': {
                'Parent': gene_id,
                'Name': f'Exon {orf_index+1}'
            }
        },
    ]
    return subfeatures

def annotate_orfs(orfs, sequence_id):
    """Annotate the identified ORFs with detailed information."""
    annotations = []
    for index, (start, stop, orf) in enumerate(orfs):
        gene_subfeatures = annotate_gene(start + 1, stop, index)
        attributes = {
            'ID': f'gene{index+1}',
            'Name': f'Gene{index+1}',
            'Alias': f'Gene{index+1}_alias',
            'Note': f'This is Gene {index+1}',
            'Ontology_term': 'GO:0003674',
            'Dbxref': f'GeneID:{index+1}',
            'Description': f'Description of Gene {index+1}',
            'Function': get_gene_function(orf),
            'Regulatory_elements': ",".join(get_regulatory_elements(orf)),
            'Features': ",".join(get_features(orf))
        }
        annotations.append({
            'Sequence_ID': sequence_id,
            'Type': 'gene',
            'Start': start + 1,
            'End': stop,
            'Strand': '+',
            'Attributes': attributes
        })
        for subfeature in gene_subfeatures:
            annotations.append({
                'Sequence_ID': sequence_id,
                'Type': subfeature['Type'],
                'Start': subfeature['Start'],
                'End': subfeature['End'],
                'Strand': '+',
                'Attributes': subfeature['Attributes']
            })
    return annotations

def save_gff3(annotations, output_file):
    """Save the annotations in GFF3 format."""
    with open(output_file, "w") as file:
        file.write("##gff-version 3\n")
        for annotation in annotations:
            attributes = ";".join([f"{key}={value}" for key, value in annotation['Attributes'].items()])
            file.write(f"{annotation['Sequence_ID']}\t.\t{annotation['Type']}\t{annotation['Start']}\t{annotation['End']}\t.\t{annotation['Strand']}\t.\t{attributes}\n")

def get_gene_function(orf):
    """Placeholder function to retrieve gene function."""
    return "Example Function"

def get_regulatory_elements(orf):
    """Placeholder function to retrieve regulatory elements."""
    return ["Regulatory Element A", "Regulatory Element B"]

def get_features(orf):
    """Placeholder function to retrieve other features."""
    return ["Feature A", "Feature B"]

def predict_genes(genome_sequence):
    """Predict genes by finding ORFs in the genome sequence."""
    gene_features = []
    orfs = find_orfs(genome_sequence)
    for start, end, _ in orfs:
        gene_location = FeatureLocation(start, end)
        gene_feature = SeqFeature(gene_location, type="gene")
        gene_features.append(gene_feature)
    return gene_features

def annotate_function(genome_features):
    """Annotate genes with sequential names."""
    for i, gene_feature in enumerate(genome_features):
        gene_feature.qualifiers["gene_name"] = f"Gene_{i+1}"
    return genome_features

def visualize_genome_annotations(genome_features, annotations, output_text):
    """Visualize annotated genome features similar to GenBank format."""
    output = ""
    output += f"LOCUS       {genome_features.id}            {len(genome_features.seq)} bp    DNA\n"
    output += f"DEFINITION  {genome_features.description}\n"
    output += f"ACCESSION   {genome_features.id}\n"
    output += "VERSION     " + genome_features.id + ".1\n"
    output += f"SOURCE      {annotations['source']}\n"
    output += f"  ORGANISM  {annotations['organism']}\n"
    output += f"            {annotations['taxonomy']}\n"
    output += "\n"
    output += "FEATURES             Location/Qualifiers\n"
    for feature in genome_features.features:
        output += f"     {feature.type}          {feature.location.start}..{feature.location.end}\n"
        for key, value in feature.qualifiers.items():
            output += f"                     /{key}=\"{value[0]}\"\n"
    output += "//"
    
    output_text.insert(tk.END, output)

def process_genome_file(genome_file, output_format, output_text):
    """Read genome file, predict and annotate genes, and visualize them."""
    genome_records = SeqIO.parse(genome_file, "fasta")
    
    for record in genome_records:
        predicted_genes = predict_genes(record.seq)
        annotated_genes = annotate_function(predicted_genes)
        annotations = {
            "source": "Gallus gallus (chicken)",
            "organism": "Gallus gallus",
            "taxonomy": "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus."
        }
        
        if output_format == "GFF3":
            annotations_gff3 = annotate_orfs(find_orfs(record.seq), record.id)
            save_gff3(annotations_gff3, "annotations.gff3")
            messagebox.showinfo("Success", "GFF3 annotations saved to annotations.gff3")
        else:
            visualize_genome_annotations(record, annotations, output_text)

def browse_genome_file():
    """Open a file dialog to select a genome file."""
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", ".fasta"), ("All files", ".*")])
    if file_path:
        output_text.delete(1.0, tk.END)
        try:
            output_format = format_var.get()
            if output_format:
                process_genome_file(file_path, output_format, output_text)
            else:
                messagebox.showerror("Error", "Please select an output format.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")

def save_output():
    """Open a file dialog to save the output."""
    file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", ".txt"), ("All files", ".*")])
    if file_path:
        with open(file_path, 'w') as file:
            file.write(output_text.get(1.0, tk.END))

# Create the main window
root = tk.Tk()
root.title("Gene Annotation Tool")

# Create and place widgets
tk.Label(root, text="Select Output Format:").pack(pady=10)
format_var = tk.StringVar(value="GFF3")
tk.Radiobutton(root, text="GFF3", variable=format_var, value="GFF3").pack()
tk.Radiobutton(root, text="GenBank", variable=format_var, value="GenBank").pack()

browse_button = tk.Button(root, text="Select Genome File", command=browse_genome_file)
browse_button.pack(pady=20)

output_text = tk.Text(root, wrap=tk.WORD, height=20, width=80)
output_text.pack(padx=20, pady=(0, 20))

save_button = tk.Button(root, text="Save Output", command=save_output)
save_button.pack(pady=10)

# Run the application
root.mainloop()
