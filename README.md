# crosslinkR

## Tasks
- read in Kojak results
- find crosslinks
- show score distribution and location of crosslinks
- calculate FDR (experimental, not used in the manuscript)

## Modified FDR calculation
A crosslink can have peptides from two different proteins. Assuming that the concentration of target and decoy hits is equal, we expect seeing target or decoy proteins half of the time each. A decoy crosslink is defined as a crosslink with one or two peptides originating from a decoy protein. Thus, the propability to see a decoy corsslink is 0.75. The calculation of the cumulative FDR needs to be corrected for this imbalance. We here propose a revised calculation:

FDR = #Decoy hits / #Target hits
FDR_corrected = (#Decoy hits / 3) / #Target hits

Compare also to the size ratio in eq. 2 from [Levitsky et al., 2016](https://pubs.acs.org/doi/pdf/10.1021/acs.jproteome.6b00144). The decoy database can be seen as thrice as large as the target database.
