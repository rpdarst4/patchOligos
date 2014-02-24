#!/usr/bin/python2.7
# (c) 2014, Russell Darst, University of Florida

"""
* retrieves patch-bisulfite oligonucleotides
* following method of Varley et al. to calculate TMs
* use TSSs=False to disable gene lookup by RefSeq ID
  (skips slow initial setup)
* use .add_locus function to add any locus to database given a unique
  identifier, chromosome, basepair, and strand
* will need to update Entrez email
* may need to change GenBankIDs (currently set for hg19)
"""

##############################################################################
##############################################################################

from Bio import Entrez,Restriction,SeqIO
from ftplib import FTP
from os import path
import gzip,sqlite3

# Entrez requires an email address for downloads
# replace with your own email address please
Entrez.email='rdarst@ufl.edu'

# names of hg19 chromosomes
GenBankIDs={'chr'+str(n): 'CM000{0}.1'.format(n+662) for n in range(1,23)}
GenBankIDs.update(chrX='CM000685.1',chrY='CM000686.1')

##############################################################################
##############################################################################

# Varley et al. method to calculate TMs
# from OligoCalc "salt adjusted",
# http://basic.northwestern.edu/biotools/OligoCalc.html
def MT(seq):
    if not seq: return -1
    w,x,y,z=[seq.upper().count(y) for y in 'ATGC']
    return 100.5+(41*(y+z)/(w+x+y+z))-(820/(w+x+y+z))-21.597097928022087

##############################################################################
##############################################################################

# creates and reads SQL file of sites and corresponding patch oligos
# so patch oligos need not be re-calculated each time
class patch_oligos():

    ##########################################################################

    # to create SQL database
    def __init__(self,db='patchOligos.db',TSSs=True):
        self.conn=sqlite3.connect(db)
        self.curs=self.conn.cursor()
        for command in (
            '{0} TABLE {1} gene (id {2}, loc, UID)',
            '{0} TABLE {1} loci (id {2}, chrom, bp, strand, versions)',
            '{0} TABLE {1} olig (id {2}, config, loc, x, y, A, B, seq)',
            '{0} TABLE {1} cfg  (id {2}, enzymes, length, TM, tries, window)',
            '{0} INDEX {1} gene_by_loc ON gene (loc)',
            '{0} INDEX {1} gene_by_UID ON gene (UID)',
            '{0} INDEX {1} loc_by_loc  ON loci (chrom, bp, strand)',
            '{0} INDEX {1} olig_by_loc ON olig (loc)'):
            self.curs.execute(command.format(
                'CREATE','IF NOT EXISTS','INTEGER PRIMARY KEY NOT NULL'))

        # default configuration
        self.curs.execute('SELECT COUNT(*) FROM cfg')
        self.__settings__=self.curs.fetchone()[0]
        if not self.__settings__:
            self.curs.execute(
                'INSERT INTO cfg VALUES (NULL, ?, ?, ?, ?, ?)',tuple([
                    str(arg) for arg in (['AluI'],(200,400),62,3,0)]))
            self.__settings__=1
        print self.settings()

        # on file creation, get all TSS positions;
        # optional, but needed for lookup by RefSeq UID
        self.curs.execute('SELECT COUNT(*) FROM gene')
        if self.curs.fetchone()[0]==0 and TSSs:
            for TSS in refGene(): self.add_locus(*TSS)
        self.conn.commit()

    ##########################################################################

    # to read SQL database (& update if needed)
    # by RefSeq UID, or any unique locus identifier
    def __getitem__(self,UID):
        
        # N = which settings
        N=self.__settings__

        # i = each row ID with UID
        for i, in self.curs.execute(
            'SELECT loc FROM gene WHERE UID = ?',(UID,)):        
            self.curs.execute(
                'SELECT chrom, bp, strand, versions FROM loci WHERE id = ?',
                (i,))
            
            # x = chromosome, y = basepair, z = strand
            x,y,z,versions=self.curs.fetchone()
            x,y,z,versions=str(x),int(y),str(z),versions.split(',')

            # if oligos have not been chosen with current settings
            if str(N) not in versions:
                if self.pick_oligos(i,x,y,z):
                    self.curs.execute(
                        'UPDATE loci SET versions = ? WHERE id = ?',
                        (','.join(versions+[str(N)]),i))
                self.conn.commit()

            # yield each option (can there be more than one?)
            for a,b,A,B,seq in self.curs.execute(
                'SELECT x, y, A, B, seq FROM olig '
                'WHERE loc = ? AND config = ?',(i,N,)):
                
                    # A, B = adaptor oligo sequences,
                    # a, b = relative positions of A and B adaptors,
                    # seq = reference sequence between adaptors
                    yield UID,x,y,z,a,b,A,B,seq

    ##########################################################################

    # store TSS positions to look-up chromosome, basepair by RefSeq UID
    def add_locus(self,UID,chrom,bp,strand):
        for i, in self.curs.execute(
            'SELECT id FROM loci WHERE chrom = ? AND bp = ? AND strand = ?',
            (chrom,bp,strand)): break
        else:
            self.curs.execute('INSERT INTO loci VALUES (NULL, ?, ?, ?, ?)',
                              (chrom,bp,strand,'0'))
            i=self.curs.lastrowid
        self.curs.execute('INSERT INTO gene VALUES (NULL, ?, ?)',(i,UID))
        self.conn.commit()

    ##########################################################################

    # method to choose oligos aka the heart of patchOligos.py
    # i = a row ID in table "loci"; x,y,z = chromosome, basepair, strand
    def pick_oligos(self,i,x,y,z):

        # read settings
        self.curs.execute(
            'SELECT * FROM cfg WHERE id = ?',(self.__settings__,))
        enzymes,length,TM,tries,window=[
            eval(item) for item in self.curs.fetchone()[1:]]
        d=max(1000,10*window)

        # try to get sequence from Entrez
        for n in range(tries):
            try: record=SeqIO.read(Entrez.efetch(
                'nucleotide',rettype='fasta',id=GenBankIDs[x],
                seq_start=y-d,seq_stop=y+d+1),'fasta')
            except: continue
            if n>1: print '{0}:{1}... {2} tries'.format(x,y,n)
            break
        else:
            print '{0}:{1}... fail'.format(x,y)
            return False
        if z=='-': record.seq=record.seq.reverse_complement()

        # store any patch-BGS adaptor candidates
        # "Restriction" module starts at 1, unlike normal Python counting
        sites=sorted(set([
            site-1 for enzyme in enzymes for site in
            Restriction.AllEnzymes.get(enzyme).search(record.seq)]))
        for a,b in zip(sites[:-1],sites[1:]):
            if not length[0]<=b-a<=length[1]: continue
            if a-window>d or b+window<d: continue
            A=B=None
            for j in range(b-a):
                if not A: A=str(record.seq[a:a+j].reverse_complement())
                if not B: B=str(record.seq[b-j:b].reverse_complement())
                if MT(A)<TM: A=None
                if MT(B)<TM: B=None
            if (A and B): self.curs.execute(
                'INSERT INTO olig VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)',
                (self.__settings__,i,a-d,b-d,A+'-adaptor','adaptor-'+B,
                 str(record.seq[a:b+1])))
        self.conn.commit()
        return True

    ##########################################################################

    # to see current settings
    def settings(self):
        self.curs.execute(
            'SELECT * FROM cfg WHERE id = ?',(self.__settings__,))
        return self.curs.fetchone()[1:]

    ##########################################################################

    # to change settings e.g. which enzyme or TM
    def update(self,**kwargs):
        update=tuple([str(val) for val in (
            kwargs.get('enzymes',['AluI']),
            kwargs.get('length',(200,400)),
            kwargs.get('TM',62),
            kwargs.get('tries',3),
            kwargs.get('window',0))])

        # check whether same settings have been used before
        for i, in self.curs.execute(
            'SELECT id FROM cfg WHERE enzymes = ? AND length = ? AND TM = ? '
            'AND tries = ? AND window = ?',update):
            self.__settings__=i
            break
        else:
            self.curs.execute(
                'INSERT INTO cfg VALUES (NULL, ?, ?, ?, ?, ?)',update)
            self.__settings__=self.curs.lastrowid
            self.conn.commit()        

##############################################################################
##############################################################################

# download TSS positions from UCSC (if not already present)
def refGene(db='/goldenPath/hg19/database/refGene.txt.gz'):
    if not path.exists(path.split(db)[1]):
        host=FTP('hgdownload.cse.ucsc.edu')
        host.login('anonymous',Entrez.email)
        host.retrbinary('RETR '+db,open(path.split(db)[1],'wb').write)
    for row in gzip.open(path.split(db)[1],'r'):
        row=row.split('\t')
        if GenBankIDs.get(row[2],False):
            yield row[1],row[2],row[4:6][row[3]=='-'],row[3]
