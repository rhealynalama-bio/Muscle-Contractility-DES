# ============================================================
# Muscle Contractility: DES 
# Transcript: NM_001927.4 | Variant: c.38C>T (p.Ser13Phe)
# ClinVar Variation ID: 44260

# ---- Step 1. Setup -----------------------------------------
# Run installation only once on your computer:
# install.packages(c("rentrez", "ggplot2"))

library(rentrez)
library(ggplot2)

gene_symbol          <- "DES"
transcript_accession <- "NM_001927.4"
clinvar_id           <- "44260"
variant_protein      <- "p.Ser13Phe"

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}

# ---- Step 2. NCBI Gene record ------------------------------
hits <- entrez_search(
  db = "gene",
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)
if (length(hits$ids) == 0) stop("Gene not found")

# Verify the exact gene symbol, not just a similar family member
gene_id <- NA_character_
for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id
    gene_info <- rec
    break
  }
}
if (is.na(gene_id)) stop("Exact human gene not found")
print(gene_info)

writeLines(c(
  paste("Gene symbol:", gene_info$name),
  paste("NCBI Gene ID:", gene_id),
  paste("Official name:", gene_info$description),
  paste("Chromosome:", gene_info$chromosome)
), "results/gene_summary.txt")

# ---- Step 3. Protein sequence from the assigned transcript -
gb <- entrez_fetch(
  db = "nuccore", id = transcript_accession,
  rettype = "gb", retmode = "text"
)

pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))
if (length(fragment) == 0) {
  stop("NP protein ID not found; inspect the transcript record")
}
protein_accession <- gsub('^/protein_id="|"$', '', fragment)

protein_fasta <- entrez_fetch(
  db = "protein", id = protein_accession,
  rettype = "fasta", retmode = "text"
)
if (!startsWith(protein_fasta, ">")) stop("Invalid FASTA")

fasta_file <- file.path("data", paste0(gene_symbol, "_protein.fasta"))
writeLines(protein_fasta, fasta_file)

cat("Matched protein:", protein_accession, "\n")
writeLines(paste("NP accession:", protein_accession),
           "results/protein_accession.txt")

# ---- Step 4. Protein length and amino-acid composition -----
fasta_lines <- readLines(fasta_file)
sequence <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))
results <- data.frame(AminoAcid = standard, Count = as.integer(counts))
write.csv(results, "results/amino_acid_composition.csv", row.names = FALSE)

cat("Protein length:", protein_length, "amino acids\n")
print(head(results[order(-results$Count), ], 3))

fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(title = paste(gene_symbol, "amino-acid composition"),
       x = "Amino acid", y = "Number of residues")
print(fig1)
ggsave("figures/01_amino_acid_composition.png", fig1,
       width = 8, height = 5, dpi = 300)

# ---- Step 6. Variant position plot (Option A) --------------
if (nzchar(variant_protein)) {
  matched  <- regmatches(variant_protein, regexpr("[0-9]+", variant_protein))
  position <- as.integer(matched)
  
  if (!is.na(position) && position >= 1 && position <= protein_length) {
    fig2 <- ggplot() +
      geom_segment(aes(x = 1, xend = protein_length, y = 0, yend = 0),
                   linewidth = 2) +
      geom_point(aes(x = position, y = 0), size = 4) +
      annotate("text", x = position, y = 0.25,
               label = variant_protein, size = 4) +
      xlim(0, protein_length) + ylim(-0.6, 0.6) +
      theme_void() +
      labs(title = paste(gene_symbol, "reference protein"),
           x = "Amino-acid position") +
      theme(axis.title.x = element_text())
    print(fig2)
    ggsave("figures/02_variant_position.png", fig2,
           width = 8, height = 2.5, dpi = 300)
  }
}

# ---- Record for README -------------------------------------
cat("\n--- Values to copy into README.md ---\n")
cat("Gene ID:", gene_id, "\n")
cat("Official name:", gene_info$description, "\n")
cat("Chromosome:", gene_info$chromosome, "\n")
cat("NP accession:", protein_accession, "\n")
cat("Protein length:", protein_length, "aa\n")
cat("Date run:", format(Sys.Date()), "\n")
sessionInfo()

# ---- Step 6b. Domain map of desmin (UniProt P17661) ----------
# Coordinates read from UniProt P17661, Family & Domains, accessed 2026-10-27
domains <- data.frame(
  Region = c("Head", "Coil 1A", "Linker 1", "Coil 1B", "Linker 12",
             "Coil 2A", "Linker 2", "Coil 2B", "Tail"),
  start  = c(2, 109, 142, 152, 253, 269, 288, 296, 413),
  end    = c(108, 141, 151, 252, 268, 287, 295, 412, 470)
)
domains$Region <- factor(domains$Region, levels = domains$Region)

variant_pos <- as.integer(regmatches(variant_protein,
                                     regexpr("[0-9]+", variant_protein)))

fig3 <- ggplot(domains) +
  geom_rect(aes(xmin = start - 1, xmax = end, ymin = 0, ymax = 1, fill = Region),
            colour = "black") +
  annotate("segment", x = variant_pos, xend = variant_pos, y = 1.6, yend = 1.05,
           arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.9) +
  annotate("text", x = variant_pos + 8, y = 1.75, hjust = 0, size = 4,
           label = paste0(variant_protein, " (residue ", variant_pos,
                          ", head domain)")) +
  scale_x_continuous(limits = c(0, 472), breaks = c(1, 100, 200, 300, 400, 470)) +
  ylim(-0.1, 2) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Human desmin (DES, P17661) domain map",
       x = "Amino-acid position (470 aa)", y = NULL, fill = "UniProt region",
       caption = "Domain coordinates: UniProt P17661, accessed 2026-10-27. Variant marker added by the student.")
print(fig3)
ggsave("figures/protein_structure_or_splice_diagram.png", fig3,
       width = 9, height = 3.5, dpi = 300)
