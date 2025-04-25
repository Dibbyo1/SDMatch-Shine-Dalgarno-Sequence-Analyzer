# SDMatch-Shine-Dalgarno-Sequence-Analyzer
SDMatch is a Python-based tool for analyzing Shine-Dalgarno sequences across upstream regions of bacterial genes. It detects mismatches, computes statistical insights, and determines the consensus sequence that likely initiates translation in Escherichia coli. </br>
This project focuses on identifying Shine-Dalgarno motifs (AGGAGG) across 20 upstream gene sequences obtained from the E. coli K12 MG1655 strain. It calculates mismatches and separation distances from the translation start site (ATG), aggregates the data in a CSV spreadsheet, and computes statistical summaries including: </br>
 - Average and standard deviation of mismatches.</br>
 - Average and standard deviation of separation distances.</br>
 - Consensus base at each position.</br>

# Features
Analyzes Shine-Dalgarno motif presence in upstream sequences.</br>
Computes per-sequence mismatch counts to the ideal AGGAGG motif.</br>
Measures separation from SD motif to translation start site.</br>
Calculates consensus motif based on most frequent base at each position.</br>
Outputs results as a clean CSV and summary statistics to terminal.</br>
