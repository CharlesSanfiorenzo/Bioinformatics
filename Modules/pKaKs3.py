#This short script uses the output values of KaKs.pl & SnpEff to calculate mutational load using Nei-Gojobori: pKa/Ks = [-3/4ln(1-4pn/3)] / [-3/4ln(1-4ps/3)], where ps = syn SNPs / syn sites and pn = nonsyn SNPs / nonsyn sites

from math import log #If for some reason you need to calculate the logarithm of a negative number, import cmath instead.
import configparser

config = configparser.RawConfigParser()
config.read("config.ini")
nonSyn_site = float(config.get("myvars", "non-synonymous_number"))
Syn_site = float(config.get("myvars", "synonymous_number"))
nonSyn_SNP = float(config.get("myvars", "non-synonymous_snp"))
Syn_SNP = float(config.get("myvars", "synonymous_snp"))

pn = nonSyn_SNP/nonSyn_site
ps = Syn_SNP/Syn_site

print("The pKs/Ks ratio for this organism is:", (-3/4*log(1-(4*pn)/3))/(-3/4*log(1-(4*ps)/3)) )
