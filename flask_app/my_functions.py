import os
from flask import jsonify, logging, make_response
import numpy as np
import pandas as pd
import dnacauldron as dc
from Bio import SeqIO
from Bio.Seq import Seq
from io import BytesIO, StringIO
from tempfile import NamedTemporaryFile, tempdir
from werkzeug.utils import secure_filename

UPLOADER_FOLDER = "/flask_app/uploads"
os.makedirs(UPLOADER_FOLDER, exist_ok=True)

tempdir = UPLOADER_FOLDER

def convert_to_csv(dfs):
    # Convert DataFrames to CSV files and store them as Python variables
    """
    Convert a list of DataFrames to CSV format.

    Args:
        dfs (list of pandas.DataFrame): List of DataFrames to be converted to CSV.

    Returns:
        list of str: List of CSV strings corresponding to the DataFrames.
    """
    csv_data = []
    for index, df in enumerate(dfs):
        csv_buffer = StringIO()
        df.to_csv(csv_buffer, index=False)
        csv_data.append(csv_buffer.getvalue())

    return csv_data

def save_uploaded_file(file):
    """
    Save an uploaded file to a temporary location.

    Args:
        file: The uploaded file object.

    Returns:
        str: The path to the saved temporary file.
    """    
    filename = secure_filename(file.filename)
    # Create a temporary file
    temp_file = NamedTemporaryFile(delete=False, prefix=filename)
    
    # Write the uploaded file to the temporary file
    file.save(temp_file.name)
    
    # Return the path of the temporary file
    return temp_file.name
############################# GET DOMESTICATED FRAGMENTS FROM GENBANK FILES ##############################################################################

def get_domesticated_fragments(domesticated_fragments_paths):
    """
    Extract domesticated fragments from GenBank files.

    Args:
        domesticated_fragments_paths (list): List of file paths to GenBank files.

    Returns:
        list: A list of dictionaries containing sequence ID, sequence, and size of each fragment.
    """
    sequences = []

    # Iterate over each file in the directory
    for file_path in domesticated_fragments_paths:
        with open(file_path, "r") as handle:
            # Parse GenBank file
            for record in SeqIO.parse(handle, "genbank"):
                # Extract sequence ID and sequence
                sequence_id = record.id
                sequence = str(record.seq)
                size = int(len(sequence)) - 2105 #remove backbone size (pUPD2)

                # Create a dictionary to hold sequence ID and sequence
                sequence_info = {"id": sequence_id, "sequence": sequence, "size": size}

                # Append dictionary to the sequences list
                sequences.append(sequence_info)
    return sequences

############################# GET BACKBONES FROM GENBANK FILES #################################################################

def get_backbones(backbones_paths):
    """
    Extract backbones from GenBank files.

    Args:
        backbones_paths (list): List of file paths to GenBank files.

    Returns:
        list: A list of GenBank records representing the backbones.
    """
    backbones = []
    # Iterate over each file in the directory
    for file_path in backbones_paths:
        with open(file_path, "r") as handle:
            # Parse GenBank file
            for record in SeqIO.parse(handle, "genbank"):
                backbones.append(record)
    return backbones

############################# GENERATING MINI-CHUNKS WITH OVERLAPS ##########################################################################

def find_and_concatenate_sequences(sequence_list, max_length, num_sequences):
    """
    Find and concatenate sequences to form larger sequences with specified constraints. This function creates overlaps between the final assembled sequences. Each overlap is one initial fragment.

    Args:
        sequence_list (list): List of dictionaries containing sequence information.
        max_length (int): Maximum allowed length for concatenated sequences.
        num_sequences (int): Number of concatenated sequences to generate.

    Returns:
        tuple: A concatenated sequence and a list of lists containing combined sequence IDs.
    """
    combined_ids_list = []
    #random.shuffle(sequence_list)
    last_sequence_id = sequence_list[0]["id"]
    last_sequence = sequence_list[0]["sequence"]
    last_sequence_length = sequence_list[0]["size"]
    
    for _ in range(num_sequences):
        concatenated_sequence2 = ""
        combined_ids2 = []
        concatenated_sequence2_size = 0

        # Add the corresponding last sequence as the first fragment of new assembly
        concatenated_sequence2 += last_sequence
        combined_ids2.append(last_sequence_id)
        concatenated_sequence2_size += last_sequence_length

        # Continue concatenating sequences until the maximum length is reached
        remaining_length = max_length - concatenated_sequence2_size
        remaining_sequences = [seq for seq in sequence_list if seq["id"] != last_sequence_id]
        #random.shuffle(remaining_sequences)

        for seq_info in remaining_sequences:
            sequence_id = seq_info["id"]
            sequence = seq_info["sequence"]
            sequence_length = seq_info["size"]

            if remaining_length >= sequence_length:
                concatenated_sequence2 += sequence
                combined_ids2.append(sequence_id)
                remaining_length -= sequence_length
                
            else:
                break
        
        combined_ids_list.append(combined_ids2)
        last_sequence_id=combined_ids2[-1]
        
        # Find the sequence corresponding to the last sequence ID
        for seq_info in sequence_list:
            if seq_info["id"] == last_sequence_id:
                last_sequence = seq_info["sequence"]
                last_sequence_length = seq_info["size"]
                break
                
        sequence_list = [seq for seq in remaining_sequences if seq["id"] not in combined_ids2] 
    
    return concatenated_sequence2, combined_ids_list


class GB_fragment:
    """
    A class to represent a Golden Braid fragment.

    Attributes:
        fragment_id (str): The ID of the fragment.
        sequence (str): The sequence of the fragment.
        backbone (str): The backbone associated with the fragment.
        size (int): The size of the fragment.
        destination (GB_fragment): The destination fragment in the assembly.
        level (numpy.array): The level of the fragment which represents the backbone. 
        origin (list): The origin fragments of the current fragment.
        order (list): The binary order of the fragment.
        decimal_order (int): The decimal order of the fragment.
        all_basal (list): The list of all basal (child) fragments.
    """
    
    def __init__(self, fragment_id="", sequence="", backbone = "", size=0, destination=None, level=np.array([0, 0]), origin=[], order=[], decimal_order = 0, all_basal = None):
        """
        Initialize a GB_fragment instance.

        Args:
            fragment_id (str): The ID of the fragment.
            sequence (str): The sequence of the fragment.
            backbone (str): The backbone associated with the fragment.
            size (int): The size of the fragment.
            destination (GB_fragment): The destination fragment in the assembly.
            level (numpy.array): The level of the fragment which represents the backbone.
            origin (list): The origin fragments of the current fragment.
            order (list): The binary order of the fragment.
            decimal_order (int): The decimal order of the fragment.
            all_basal (list): The list of all basal (child) fragments.
        """
        self.fragment_id = fragment_id
        self.sequence = sequence
        self.backbone = backbone
        self.size = size
        self.level = level
        self.origin = origin
        self.destination = destination
        self.order = order
        self.all_basal = [] #children
        self.decimal_order = decimal_order
        
    def start_forwardpropagation(self):
        """Initiate forward propagation to update destination fragment."""
        self.destination.forwardpropagation()

    def forwardpropagation(self):
        """Propagate sequence and size information from origin fragments."""
        self.sequence = self.origin[0].sequence + self.origin[1].sequence
        self.size = self.origin[0].size + self.origin[1].size
        self.all_basal = self.origin[0].all_basal + self.origin[1].all_basal
        if self.destination:
            self.destination.forwardpropagation()

    def backpropagation(self):
        """Propagate level and order information back to origin fragments."""
        if self.origin:
            changes = [(1,0), (1,1)]
            changes2 = [(1,1), (1,0)]
            for i in range(2):
                if self.level[1] == 0:
                    self.origin[i].level = np.mod(self.level + changes[i], 2)
                else:
                    self.origin[i].level = np.mod(self.level + changes2[i], 2)
                self.origin[i].order = self.order.copy()
                self.origin[i].order.append(self.origin[i].level[1])
                self.origin[i].backpropagation()
    
    def define_backbone(self):
        """
        Define the backbone for the fragment based on its level.

        Returns:
            numpy.array: The level of the fragment.
        """
        if np.array_equal(self.level, np.array([0,0])):
            self.backbone = "pdgb3_alpha1"
        elif np.array_equal(self.level, np.array([0,1])):
            self.backbone = "pdgb3_alpha2"
        elif np.array_equal(self.level, np.array([1,0])):
            self.backbone = "pdgb3_omega1"
        elif np.array_equal(self.level, np.array([1,1])):
            self.backbone = "pdgb3_omega2"
            
        return self.level
    
############################# INITIATION OF ASSEMBLY TREE ###################################################################################

def assemble (plasmid1, plasmid2):
    """
    Assemble two plasmids into a new GB_fragment.

    Args:
        plasmid1 (GB_fragment): The first plasmid fragment.
        plasmid2 (GB_fragment): The second plasmid fragment.

    Returns:
        GB_fragment: The assembled plasmid fragment.
    """
    plasmid3 = GB_fragment(origin = [plasmid1, plasmid2], all_basal= plasmid1.all_basal + plasmid2.all_basal)
    plasmid1.destination = plasmid3
    plasmid2.destination = plasmid3
    return plasmid3
    
def convert_binary(initial, initial_set):
    """
    Convert the binary order of a fragment to a decimal value.

    Args:
        initial (GB_fragment): The fragment to convert.
        initial_set (list): The set of initial fragments.

    Returns:
        int: The decimal order of the fragment.
    """
    generation = max(len(fragment.order) for fragment in initial_set)
    decimal_value = 0
    adj_order = initial.order + [0]*(generation - len(initial.order))
    for gen in range(0, generation):
        decimal_value = (decimal_value << 1) | adj_order[gen]
    return int(decimal_value)

def calculate_alpha_probability(n):
    """
    Calculate the probability to have an alpha backbone for a given number of fragments.

    Args:
        n (int): The number of fragments.

    Returns:
        float: The alpha probability.
    """
    alpha_probability = (0.5*n+(1.5*n-2**np.ceil(np.log2(n)))*(-1)**np.mod(np.ceil(np.log2(n)), 2))/n
    return alpha_probability


def build_tree(combined_ids_list, sequences):
    """
    Build an assembly tree from combined IDs and sequences.

    Args:
        combined_ids_list (list): List of combined sequence IDs.
        sequences (list): List of sequence information dictionaries.

    Returns:
        list: A list of DataFrames representing the assembly plan for each final sequence.
    """
    dfs = []
    for p in range (0, len(combined_ids_list)):
        names = combined_ids_list[p]
        fragments = []
        fake_fragments = []

        for name in names:
            fake = GB_fragment(fragment_id = name)
            fake.all_basal.append(fake.fragment_id)
            fake_fragments.append(fake)
            for seq_info in sequences:
                sequencee = seq_info['sequence']
                sequence_id = seq_info['id']
                sequence_size = seq_info["size"]
                if name == sequence_id:
                    fragments.append(seq_info)
                    break

        initial_set = fake_fragments.copy()
        merged_fragments = []
        all_assembled = []
        all_orders = []
        all_decimal_orders = []

        while len(fake_fragments) > 1:
            if len(fake_fragments) % 2 != 0:
                last_fragment = fake_fragments.pop()
                merged_fragments.insert(0, last_fragment)
    
            for i in range (0, len(fake_fragments), 2):
                plasmid1 = fake_fragments[i]
                plasmid2 = fake_fragments[i + 1]
                plasmid3 = assemble(plasmid1, plasmid2) 
        
                merged_fragments.append(plasmid3)
                all_assembled.append(plasmid3)

            fake_fragments = merged_fragments.copy()
            merged_fragments = []  

    
        alpha_probability = calculate_alpha_probability(len(initial_set))

        fragment_counter = 0

        for fragment in fake_fragments:
            if alpha_probability < 0.5:
                fragment.level = np.array([1, 0])
            fragment.backpropagation()

        for initial in initial_set:
            initial.decimal_order = convert_binary(initial, initial_set)
            all_orders.append(initial.decimal_order)

        sorted_ordr = sorted(all_orders)
        # Adjust the indices to maintain the initial order
        adjusted_ordr = [sorted_ordr.index(value) for value in all_orders]
        
        for initial in initial_set:
            initial.decimal_order = adjusted_ordr.pop(0)
            all_decimal_orders.append(initial.decimal_order)

        reordered_fragments = [fragments[v] for v in all_decimal_orders]

        for initial, seq_info in zip(initial_set, reordered_fragments):
            sequencee = seq_info['sequence']
            sequence_id = seq_info['id']
            sequence_size = seq_info["size"]
            initial.fragment_id = sequence_id
            initial.sequence = sequencee
            initial.size = sequence_size
                    
            initial.define_backbone()
            initial.start_forwardpropagation()

        for assembled in all_assembled:
            assembled.define_backbone()
            fragment_counter += 1
            assembled.fragment_id =  f"construct_{fragment_counter}"

        max_parts = 4

        # Initialize the DataFrame columns
        columns = ['construct'] + ['parts'] + [''] + [' ']

        # Initialize an empty list to store row data
        data = []

        # Populate the row data
        for initial in initial_set:
            row = [initial.fragment_id.split('_')[1]] + [initial.fragment_id] + [initial.backbone]
            # Pad with None if the row has fewer parts than max_parts
            row += [None] * (max_parts - len(row))
            data.append(row)
                
        for assembled in all_assembled:
            row = [split_if(assembled.fragment_id)] + [split_if(assembled.origin[0].fragment_id)] + [split_if(assembled.origin[1].fragment_id)] + [assembled.backbone]
            
            # Pad with None if the row has fewer parts than max_parts
            row += [None] * (max_parts - len(row))
            data.append(row)

        # Create the DataFrame
        df = pd.DataFrame(data, columns=columns)
        dfs.append(df)

    return dfs

def split_if (fragment_id):
    """
    Split the fragment ID if it contains 'pUU'.

    Args:
        fragment_id (str): The fragment ID.

    Returns:
        str: The split fragment ID or the original fragment ID.
    """
    if "pUU" in fragment_id:
        return fragment_id.split('_')[1]
    return fragment_id

############################ PRIMER DESIGN ############################
import os
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import csv

def calculate_tm(primer):
    """Calculate the melting temperature of the primer."""
    return mt.Tm_NN(primer)

def gc_content(primer):
    """Calculate the GC content of the primer."""
    g = primer.count('G')
    c = primer.count('C')
    return (g + c) / len(primer) * 100

def find_primers(sequence, start, end, optimal_tm=58, tm_range=(55, 60), primer_length_range=(15, 20), gc_content_range=(40, 60)):
    """Find primers within the specified range and conditions."""
    best_primer = None
    best_tm_diff = float('inf')

    for length in range(primer_length_range[0], primer_length_range[1] + 1):
        for i in range(start, end - length + 1):
            primer = sequence[i:i + length]
            if primer[-2:] in ['GC', 'CG']:
                tm = calculate_tm(primer)
                gc = gc_content(primer)
                if tm_range[0] <= tm <= tm_range[1] and gc_content_range[0] <= gc <= gc_content_range[1]:
                    tm_diff = abs(optimal_tm - tm)
                    if tm_diff < best_tm_diff:
                        best_tm_diff = tm_diff
                        best_primer = primer

    return best_primer

def design_primers(sequence, optimal_tm=58):
    """Design forward and reverse primers for the given sequence."""
    seq_len = len(sequence)
    forward_primer = find_primers(sequence, 950, 1300, optimal_tm)
    reverse_primer = find_primers(sequence.reverse_complement(), 950, 1300, optimal_tm)
    
    if forward_primer and reverse_primer:
        return forward_primer, reverse_primer
    
    return None, None

def process_genbank_files(gb_dir):
    """Process GenBank files and design primers for each sequence."""
    primer_designs = []

    for filename in os.listdir(gb_dir):
        if filename.endswith('.gb') or filename.endswith('.gbk'):
            file_path = os.path.join(gb_dir, filename)
            record = SeqIO.read(file_path, "genbank")
            sequence = record.seq
            forward_primer, reverse_primer = design_primers(sequence)
            
            if forward_primer and reverse_primer:
                primer_designs.append({
                    'id': record.id,
                    'forward_primer': str(forward_primer),
                    'reverse_primer': str(reverse_primer)
                })

    return primer_designs

def save_to_csv(primer_designs, output_file='primer_designs_outwards.csv'):
    """Save the primer designs to a CSV file."""
    with open(output_file, mode='w', newline='') as csv_file:
        fieldnames = ['id', 'forward_primer', 'reverse_primer']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        
        writer.writeheader()
        for design in primer_designs:
            writer.writerow(design)


def generate_reports(csv_data, domesticated_paths, backbones_paths):
    """
    Generate assembly reports from CSV data.

    Args:
        csv_data (list of str): List of CSV strings representing assembly data.
        domesticated_paths (list of str): Paths to domesticated fragment files.
        backbones_paths (list of str): Paths to backbone files.

    Returns:
        flask.Response: A Flask response object with the report file for download, or None if an error occurs.

    Raises:
        Exception: If an error occurs during the report generation process.
    """

    report_paths = []
    report_dir = './shared/output/assembly_reports'
    os.makedirs(report_dir, exist_ok=True)

    for index, data in enumerate(csv_data):
        if not data:
            #g.info(f"Skipping assembly simulation for DataFrame {index + 1} due to empty data")
            continue

        #logging.info(f"Running assembly simulation for DataFrame {index + 1}")

        df = pd.read_csv(StringIO(data))
        repository = dc.SequenceRepository()
        repository.import_records(files=domesticated_paths, use_file_names_as_ids=True)
        repository.import_records(files=backbones_paths, use_file_names_as_ids=True)

        omegas_list = []
        level_0s = df.iloc[:,-1].isna().sum()

        for indexx, value in df.iloc[:level_0s, 2].items():
            if 'omega' in value.lower():
                omegas_list.append(df.at[indexx, 'parts'])
        #This if performs a change of recognision site as a quick fix regarding the lack of support for multiple restriction enzymes in DNA cauldron
        #In this way we use a single enzyme, the user does not see this hack
        if omegas_list:
            for x in omegas_list:
                sequence = str(repository.get_record(x).seq)
                sequence = sequence.replace('GGTCTC', 'CGTCTC')
                sequence = sequence.replace('GAGACC', 'GAGACG')
                modified_sequence = Seq(sequence)
                repository.get_record(x).seq = modified_sequence

        assembly_plan = dc.AssemblyPlan.from_spreadsheet(
            dataframe=df,
            is_csv=None,
            sheet_name=None,
            assembly_class=dc.Type2sRestrictionAssembly
        )

        plan_simulation = assembly_plan.simulate(sequence_repository=repository)
        report_path = os.path.join(report_dir, f"assembly_report_{index}")
        report_writer = dc.AssemblyReportWriter(
            include_mix_graphs=True, include_assembly_plots=True
        )
        plan_simulation.write_report(report_path, assembly_report_writer=report_writer)
        #send_file(report_path, as_attachment=True, download_name=os.path.basename(report_path))
        #report_paths.append(report_path)

        #shutil.copy(report_path, f'/app/assembly_reports/{os.path.basename(report_path)}')

        #for report_path in report_paths:
            #return send_file(report_path, as_attachment=True, download_name=os.path.basename(report_path))

    return report_paths
    

def run_assembly_simulation(csv_data, domesticated_paths, backbones_paths, index):
    """
    Run the assembly simulation and generate reports.

    Args:
        csv_data (list of str): List of CSV strings representing assembly data.
        domesticated_paths (list of str): Paths to domesticated fragment files.
        backbones_paths (list of str): Paths to backbone files.

    Raises:
        Exception: If an error occurs during the assembly simulation process.
    """

    report_paths = generate_reports(csv_data, domesticated_paths, backbones_paths, index)

    return report_paths # Store report paths in assembly_status

