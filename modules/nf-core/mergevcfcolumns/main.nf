process MERGE_VCF_COLUMNS {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"
    
    input:
    tuple val(meta), path(merged_tab)  // Concatenated tab file from previous step
    tuple val(meta), path(vcf), path(tbi)                           // Raw VCF file

    output:
    tuple val(meta), path("${meta.id}_final.tab.gz")   , emit: cleaned_output

    script:
    """
    python3 <<CODE
    import gzip
    import pandas as pd

    # Function to normalize 'Location' by extracting the first position (pos1) if allele is '-', 
    # subtracting 1, and returning 'chrom:new_pos' format.
    # 'chrom:new_pos' format will be used to match VCF
    def normalize_annotation_location(loc, allele):

        # Normalize 'Location' by extracting 'chrom:pos1', subtracting 1 if 'Allele' is '-', and returning the adjusted 'chrom:new_pos' format.
        
        # Parameters:
        #    loc (str): Location string in the format 'chrom:pos1-pos2'.
        #    allele (str): Allele string to determine if subtraction is needed.
        
        # Returns:
        #    str: Adjusted 'chrom:pos' string.

        if not isinstance(loc, str):
            print(f"Unexpected value type: {type(loc)} - {loc}")
            return loc  # Return original if it's not a string
        
        # Extract "chrom:pos1" before '-'
        chrom_pos = loc.split('-')[0]
        
        if ':' in chrom_pos:
            chrom, pos1 = chrom_pos.split(':')
            try:
                pos1 = int(pos1)  # Convert position to integer
                
                # Subtract 1 from position **only if Allele is '-'**
                if allele == '-':
                    pos1 -= 1
                
                return f"{chrom}:{pos1}"
            except ValueError:
                print(f"Warning: Invalid position format in Location: {loc}")
                return loc
        else:
            print(f"Warning: Unexpected format in Location: {loc}")
            return loc

    # Define the file paths
    vcf_file_path = '${vcf}'
    annotated_file_path = '${merged_tab}'
    output_file_path = '${meta.id}_merged_with_vcf_data.tab.gz'

    # Load the VCF file and dynamically get the 10th column name
    vcf_columns = []
    vcf_data = []
    with gzip.open(vcf_file_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                vcf_columns = line.strip().split('\t')
                sample_names = vcf_columns[9:]  # Get all sample names (trio support)
                break

    # Now read the VCF data lines and extract relevant columns
    with gzip.open(vcf_file_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            fields = line.strip().split('\\t')

            # ✅ Define vcf_entry before using it
            vcf_entry = {
                '#CHROM': fields[0],
                'POS': fields[1],
                'QUAL': fields[5],
                'INFO': fields[7],
                'FORMAT': fields[8]
            }

            # Add all trio sample data
            for idx, sample in enumerate(sample_names):
                vcf_entry[sample] = fields[9 + idx]
            
            vcf_data.append(vcf_entry)

    # Convert to DataFrame
    vcf_df = pd.DataFrame(vcf_data, columns=['#CHROM', 'POS', 'QUAL', 'INFO', 'FORMAT'] + sample_names)
    # Create Location column in format '#CHROM:POS'
    vcf_df['Location'] = vcf_df['#CHROM'] + ':' + vcf_df['POS']

    # Find and load the annotated file header row
    with gzip.open(annotated_file_path, 'rt') as f:
        header_row = 0
        for line in f:
            if line.startswith('#') and not line.startswith('##'):
                header = line.strip().split('\\t')  # Keep this header line
                break
            header_row += 1

    # Now load the data from the annotated file, skipping all lines up to the found header
    annotated_df = pd.read_csv(
        annotated_file_path,
        sep='\\t',
        compression='gzip',
        header=header_row,
        low_memory=False
    )
    annotated_df.columns = header  # Set the header correctly

    # Check if 'Location' exists in annotated_df for merging
    if 'Location' not in annotated_df.columns:
        print("The 'Location' column is not found in annotated_df. Please check the column names.")
        exit(1)

    # ✅ Apply the function correctly by passing both 'Location' and 'Allele'
    annotated_df['Normalized_Location'] = annotated_df.apply(
        lambda row: normalize_annotation_location(row['Location'], row['Allele']), axis=1
    )

    # Merge data based on the Location column
    merged_df = pd.merge(annotated_df, vcf_df[['Location', 'QUAL', 'INFO', 'FORMAT'] + sample_names], 
                    left_on='Normalized_Location', right_on='Location', how='left')

    # Fill missing values with "./.:.:.:.:."
    for sample in sample_names:
        merged_df[sample] = merged_df.apply(
            lambda row: row[sample] if row['IND'] == sample else "./.:.:.:.:.",
            axis=1
        )

    # Replace missing values with "-"
    merged_df = merged_df.fillna("-")

    # Remove unnecessary duplicate location and rename column
    merged_df = merged_df.drop(columns=['Normalized_Location']).rename(columns={'Location_y': 'Chrom_Pos', 'Location_x': 'Location'})

    # Determine if trio case
    if len(sample_names) == 3:
        # Select the value that matches 'IND' column
        merged_df["Sample"] = merged_df.apply(
            lambda row: row[row["IND"]] if row["IND"] in sample_names and row[row["IND"]] != "./.:.:.:.:." else "-", 
            axis=1
        )
        # Drop the original trio columns
        merged_df = merged_df.drop(columns=sample_names)

    # Save the result to a new file
    merged_df.to_csv(output_file_path, sep='\t', index=False, compression='gzip')
    print("Merged file saved to:", output_file_path)
    matched_rows = merged_df[merged_df['QUAL'].notna()].shape[0]
    total_rows = merged_df.shape[0]
    print(f"Merging complete: {matched_rows}/{total_rows} rows matched.")
    CODE

    # Run zcat, grep, awk, and bgzip to clean up duplicates in the merged file
    zcat ${meta.id}_merged_with_vcf_data.tab.gz | awk '!/^##/ && (!/^#/ || !seen++)' | bgzip > ${meta.id}_final.tab.gz
    """
}
