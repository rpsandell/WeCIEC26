# Data used in the paper can be obtained by running the following code in the R statistical computing environment.

text<-scan(file="LIVRootsAll.txt", what="char", sep="\n") # Read in the full file of 1181 roots from the LIV's reverse index.
text.sub <- gsub("(.?.?.? )?(.?.?.? )?(.?.?.? )?(.) ([ea]) (.)( .?.?.?)?( .?.?.?)?( .?.?.?)?", "\\4 \\5 \\6", text) # This regex mostly correctly captures all and only CVC sequences
CeC.table <- table(text.sub) # Reduce the resulting vector to a table of types with their frequencies
CeC.frame <- as.data.frame(CeC.table) # convert that table to a dataframe
write.table(CeC.frame, "CeCPatterns.txt", row.names=F) # write out that dataframe to a .txt file
# Editing of the output here by hand: open the file CeCPatterns.txt with any plain text editor
  # The following roots are deleted because they are not unique types:
  # The following outputs were adjusted to reduce them to the desired /CVC/ shape:
  # The following forms may or may not represent a unique type, and therefore were counted as less than 1 in proportion to the raw likelihood that they represent a unique type.


#Working on just roots of the shape /CVC-/
bicons.roots <- grep("^.[cjwh]?[wh]? [ea] .[cjwh]?[wh]?$", text) # regex finds all biconsonantal roots (^ and $ indicate beginning and end of line)
bicons.roots.results <- text[bicons.roots]
bicons.table <- table(bicons.roots.results) # reduce to table of types
bicons.frame <- as.data.frame(bicons.table) # data frame contains 182 items. Some surely non-unique, estimate at 180, since that accords with an earlier hand count that I made.
bicons.matrix <- matrix(c(1,25,180,625), nrow=2)
chisq.test(bicons.matrix)

# Working on all -CVC- sequences within roots

CeC.frame2 <- read.table("CeCPatterns.txt", sep=" ") # read back in the edited text file
# The data frame, ignoring the first row, has 308 items
.3+.5+.5+.2+.6+.5+.3+.5+.5+.5+.3+.3+.5+.3+.5 # The values for possibly unique items (= 6.3, round down to 6)

CeC.matrix <- matrix(c(4, 25, 314, 625), nrow=2) # make a matrix for CiVCi items. 4 types of /CiVCi/, 25 possible; 314 types of /CVC/, 325 possible
chisq.test(CeC.matrix) # run chi-squared test

CeC.shapes <- CeC.frame2$V1 # extract a vector of all of the /CVC/ shapes from the CeC.frame2
CeC.shapes <- CeC.shapes[2:324] # eliminate the unwanted item at the top

TeR.all <- grep("[^mnlrywsHhxqw][wh]?[wh]? [ea] [mnlryw]", CeC.shapes) # regex finds all patterns of the form /TVR/
TeR.all.results <- CeC.shapes[TeR.all]
TeR.results.frame <- as.data.frame(TeR.results.all)
write.table(TeR.results.frame, "TeRresults.txt", row.names=FALSE)
15*6 # = 90, the number of possible fixed order combinations of a stop (t, d, dh, p, b, bh, k, g, gh, c, j, jh, kw, gw, gwh) and a sonorant (m, n, r, l, y, w)
TeR.matrix <- matrix(c(77,90,314,625), nrow=2)
chisq.test(TeR.matrix) # run chi-squared test


# Not included in the paper, but sequences of two obstruents are also underattested /TeT/, two stop sequences. One might also consider adding here any TT sequences (e.g., /dh gwh-/) that occur.
TeT.all <- grep("[^mnlrywRNsHhxqw][wh]?[wh]? [ea] [^mnlrywRNsHhxqw][wh]?[wh]?", CeC.shapes)
TeT.all.results <- CeC.shapes[TeT.all]
15^2 # = 225, the number of possible permutations of the 16 stop consonants (t, d, dh, p, b, bh, k, g, gh, c, j, jh, kw, gw, gwh), allowing for replacement
TeT.frame <- as.data.frame(TeT.all.results) # 41 items --  stop combinations in general are almost as badly underattested as CiVCi!
write.table(TeR.results.frame, "TeRresults.txt", row.names=FALSE)
TeT.matrix <- matrix(c(41, 225, 314, 625), nrow=2)
chisq.test(TeT.matrix) # run chi-squared test. X-squared = 31.5499, df = 1, p-value = 1.944e-08

# Looking at roots with the same place and manner of articulation
labial <- grep("[pb]h? [ea] [pb]h?", CeC.shapes) # labial stops
labial.results <- CeC.shapes[labial] # NONE; 9 possible permutations
dental <- grep("[td]h? [ea] [td]h?", CeC.shapes) # dental stops
dental.results <- CeC.shapes[dental] # "t e t"; 9 possible permutations
palatal <- grep("[(kc)cj(gj)G]h? [ea] [(kc)cj(gj)G]h?", CeC.shapes) # palatovelar stops
palatal.results <- CeC.shapes[palatal] # "k e k"   "kc a gh"; 9 possible permutations
velar <- grep("[(kc)kg(gj)KG]w?h? [ea] [(kc)kg(gj)KG]w?h?", CeC.shapes)
velar.results <- CeC.shapes[velar] # returns "c e Kw"  "Gw e jh" "k e k"   "kc a gh" "kw e c", but visual inspection lets us eliminate 3 of these.
                                    # 36 possible permutations
# m, n, s, w, y are all segments with unique place and manner of articulation. These were checked in previous /CVC/ evaluation. Only /ses/ is found.
liquids <- grep("[Rrl] [ea] [Rrl]", CeC.shapes) 
liquids.results <- CeC.shapes[liquids] # None # 4 permutations
laryngeals <- grep("^[Hhqx] [ea] [Hhqx]$", CeC.shapes) # Need beginning ^ and ending $ to avoid catching aspirated stops
laryngeals.results <- CeC.shapes[laryngeals] # "h e h" "x e h" "x e q"; 9 permutations
    # All permutations sum to 9+9+9+36+1+1+1+1+1+4+9 = 81
plc.man.matrix <- matrix(c(6.5, 81, 314, 625), nrow=2)
chisq.test(plc.man.matrix)

# CeRC and CReC for shared place and manner
CReC.labial <- grep("[pb]h? [mnRlrwy] [ea] [pb]h?", text)
CReC.labial.results <- text[CReC.labial] # "p r e p" 9*6 permutations
CReC.dental <- grep("[td]h? [mnRlrwy] [ea] [td]h?", text)
CReC.dental.results <- text[CReC.dental] # None 9*6 permutations
CReC.palatal <- grep("[(kc)cj(gj)G]]h? [mnRlrwy] [ea] [(kc)cj(gj)G]]h?", text)
CReC.palatal.results <- text[CReC.palatal] # None 9*6 permutations
CReC.velar <- grep("[(kc)kg(gj)KG]w?h? [mnRlrwy] [ea] [(kc)kg(gj)KG]w?h?", text)
CReC.velar.results <- text[CReC.velar] # None, *6 permutations
# It is easily checked that where C = a sonorant, no CReC sequence exist EXCEPT for /wrey-/ and /wley-/
          # m*6 =6, n*6=6, r/l 4*6=24, w*6=6, y*6=6
CReC.laryngeals <- grep("^[Hhqx] [mnRNlrwy] [ea] [Hhqx]$", text)
laryngeals.results <- text[CReC.laryngeals] # "H y e h" "x l e h" "x m e h" "h r e h" "x w e h" "q n e x" "h w e x" "x n e q" 
"h w e x" 9*6. Surpisingly many.
                                            # If h1 has a different place feature, then there are probably only 3 here. 

CeRC.labial <- CReC.labial <- grep("[pb]h? [ea] [mnRlrwy] [pb]h?", text)
CeRC.labial.results <- text[CeRC.labial] # None
CeRC.dental <- grep("[td]h? [ea] [mnRlrwy] [td]h?", text)
CeRC.dental.results <- text[CReC.dental] # "t e n d"   "t e r d"   "s t e w d" "d e r dh"  "dh e w dh"
CeRC.palatal <- grep("[(kc)cj(gj)G]]h? [ea] [mnRlrwy] [(kc)cj(gj)G]]h?", text)
CeRC.palatal.results <- text[CeRC.palatal] # None
CeRC.velar <- grep("[(kc)kg(gj)KG]w?h? [ea] [mnRlrwy] [(kc)kg(gj)KG]w?h?", text)
CeRC.velar.results <- text[CeRC.velar] # "s k e n gj"  "gjh e y gjh" "Kw e w k" 
CeRC.laryngeals <- grep("^[Hhqx] [ea] [mnRNlrwy] [Hhqx]$", text)
CeRC.laryngeals.results <- text[CeRC.laryngeals] # "x e m H" "h e r H" "x e r H" "h e w H" "x e w H" "q e l h" "x e n h" "h e l x" "x e l x" "h e r x" "h e w x" "x e m q" "x e r q"

# Possible number of CReC and CeRC where C is same place and manner is 81*12; total possible CeRC and CReC permutations is 625*12

# Finding all CRVC and CVRC roots
CReC.roots <- grep("(..?.? )([mnRNlrwy] )([ea])( ..?.?)", text)
CReC.roots.results <- text[CReC.roots]
CReC.roots.reduced <- gsub("(.*)(\\S\\S?\\S?) ([mnRNlrwy]) ([ea]) (\\S\\S?\\S?)(.*)", "\\2 \\3 \\4 \\5", CReC.roots.results)
CeRC.roots <- grep("(..?.? )([ea])( [mnRNlrwy])( ..?.?)", text)
CeRC.roots.results <- text[CeRC.roots]
CeRC.roots.reduced <- gsub("(.*)(\\S\\S?\\S?) ([ea]) ([mnRNlrwy]) (\\S\\S?\\S?)(.*)", "\\2 \\3 \\4 \\5", CeRC.roots.results)

CReC.roots.table <- table(CReC.roots.reduced)
CReC.roots.frame <- as.data.frame(CReC.roots.table) #226 rows
CeRC.roots.table <- table(CeRC.roots.reduced)
CeRC.roots.frame <- as.data.frame(CeRC.roots.table) # 426 rows

226+426-15 # = 637 15 estimate of number of false novelties caused by uncertain reconstructions

CRC.matrix <- matrix(c(32, 972, 637, 7500), nrow=2)
chisq.test(CRC.matrix)