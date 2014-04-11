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
import gzip,sqlite3,sys

# Entrez requires an email address for downloads
# replace with your own email address please
Entrez.email='rdarst@ufl.edu'

# Support for multiple organisms (added by A.Riva)

GenBankIDs = {}
RefGenePaths = {'hs': '/goldenPath/hg19/database/refGene.txt.gz',
                'mm': '/goldenPath/mm10/database/refGene.txt.gz',
                'mm9': '/goldenPath/mm9/database/refGene.txt.gz'}

def initGenBankIDs(org="hs"):
    global GenBankIDs

    if org == "hs":
        # names of hg19 chromosomes
        GenBankIDs = {'chr'+str(n): 'CM000{0}.1'.format(n+662) for n in range(1,23)}
        GenBankIDs.update(chrX='CM000685.1',chrY='CM000686.1')
    elif org == "mm":
        # names of mm10 chromosomes
        GenBankIDs = {'chr'+str(n): 'CM00{0}.2'.format(n+993) for n in range(1,20)}
        GenBankIDs.update(chrX='CM001013.2',chrY='CM001014.2')
    elif org == "mm9":
        # names of mm9 chromosomes
        GenBankIDs = {'chr'+str(n): 'CM00{0}.1'.format(n+993) for n in range(1,20)}
        GenBankIDs.update(chrX='CM001013.1',chrY='CM001014.1')
        
        

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

    organism = "hs"             # Keep track of which organism we're working with
    deferCommit = False;        # To speed up parsing of refseq table

    ##########################################################################

    # to create SQL database
    def __init__(self,org='hs',db='patchOligos.db',TSSs=True):
        self.organism = org
        initGenBankIDs(org=org)
        db = org + "-" + db
        self.conn=sqlite3.connect(db)
        self.curs=self.conn.cursor()
        for command in (
            '{0} TABLE {1} gene (id {2}, loc, UID)',
            '{0} TABLE {1} loci (id {2}, chrom, bp, strand, versions)',
            '{0} TABLE {1} olig (id {2}, config, loc, x, y, A, B, seq, MTA, MTB)',
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
        # print self.settings()

        # on file creation, get all TSS positions;
        # optional, but needed for lookup by RefSeq UID
        self.curs.execute('SELECT COUNT(*) FROM gene')
        if self.curs.fetchone()[0]==0 and TSSs:
            self.deferCommit = True
            for TSS in refGene(org): self.add_locus(*TSS)
            self.deferCommit = False
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
            # print("N={}, versions={}".format(N, versions))
            if str(N) not in versions:
                if self.pick_oligos(i,x,y,z):
                    self.curs.execute(
                        'UPDATE loci SET versions = ? WHERE id = ?',
                        (','.join(versions+[str(N)]),i))
                self.conn.commit()

            # yield each option (can there be more than one?)
            for a,b,A,B,seq,MTA,MTB in self.curs.execute(
                'SELECT x, y, A, B, seq, MTA, MTB FROM olig '
                'WHERE loc = ? AND config = ?',(i,N,)):
                
                    # A, B = adaptor oligo sequences,
                    # a, b = relative positions of A and B adaptors,
                    # seq = reference sequence between adaptors
                    # MTA, MTB = A and B melting temperatures
                    yield UID,x,y,z,a,b,A,B,seq,MTA,MTB

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
        if not self.deferCommit:
            self.conn.commit()

    ##########################################################################

    # method to choose oligos aka the heart of patchOligos.py
    # i = a row ID in table "loci"; x,y,z = chromosome, basepair, strand
    def pick_oligos(self,i,x,y,z):

        # print("picking oligos for {}...".format(i))
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
            # print record
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
            MTA=MTB=0
            for j in range(b-a):
                if not A: A=str(record.seq[a:a+j].reverse_complement())
                if not B: B=str(record.seq[b-j:b].reverse_complement())
                MTA = MT(A)
                MTB = MT(B)
                if MT(A)<TM: 
                    A=None
                    MTA=0
                if MT(B)<TM: 
                    B=None
                    MTB=0
            if (A and B): self.curs.execute(
                    'INSERT INTO olig VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                    (self.__settings__,i,a-d,b-d,A+'-adaptor','adaptor-'+B, str(record.seq[a:b+1]), MTA, MTB))
        self.conn.commit()
        return True

    ##########################################################################

    # to retrieve current settings
    # if the first argument is specified, select the corresponding settings
    def settings(self, idx=False):
        if idx:
            self.__settings__=idx
        else:
            idx=self.__settings__
        self.curs.execute(
            'SELECT * FROM cfg WHERE id = ?',(idx,))
        return self.curs.fetchone()[1:]

    # to show current settings
    def showSettings(self):
        st = self.settings()
        print("Enzymes: " + st[0])
        print("Region:  " + st[1])
        print("TM:      " + st[2])
        print("Tries:   " + st[3])
        print("Window:  " + st[4])

    # to list all saved settings
    def listSettings(self):
        for row in self.curs.execute("SELECT * FROM cfg ORDER BY id"):
            print row

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
def refGene(org='hs'):
    db = RefGenePaths[org]
    remoteName = path.split(db)[1]
    localPath = org + "-" + remoteName
    if not path.exists(localPath):
        print("Downloading {} to {}...".format(db, localPath))
        host = FTP('hgdownload.cse.ucsc.edu')
        host.login('anonymous', Entrez.email)
        host.retrbinary('RETR ' + db, open(localPath, 'wb').write)
    print("Parsing {}...".format(localPath))
    for row in gzip.open(localPath,'r'):
        row = row.split('\t')
        if GenBankIDs.get(row[2],False):
            yield row[1],row[2],row[4:6][row[3]=='-'],row[3]

##############################################################################
# Following functions added by A.Riva 
##############################################################################

def processFile(infile, outfile, column=0, org='hs', cfg=False):
    DB = patch_oligos(org)
    DB.settings(idx=cfg)
    with open(outfile,'w') as handle:
        handle.write('UID\tchromosome\tbasepair\tstrand\tA\tB\tsequence\tMT(A)\tMT(B)\n')
        for line in open(infile,'r'):
            gene = line.rstrip('\n').split('\t')[column]
            print("Designing primers for {}...".format(gene))
            for i in DB[gene]:
                handle.write('\n'+'\t'.join([str(j) for j in i]))

# Main function (calls processFile)
if __name__ == "__main__":

    org = "hs"
    column = 0
    sidx = False                # Use default settings
    infile = False
    outfile = False

    # Parse command-line arguments
    next = False
    for arg in sys.argv[1:]:
        if arg == "-o": next = "o"
        elif arg == "-c": next = "c"
        elif arg == "-s": next = "s"
        else:
            if next == "o": org = arg
            elif next == "c": column = int(arg)
            elif next == "s": sidx = int(arg)
            elif infile == False:
                infile = arg
            else:
                outfile = arg
            next = False

    if infile and outfile:
        processFile(sys.argv[1], sys.argv[2], org=org, column=column, cfg=sidx)
    else:
        DB = patch_oligos(org)
        print """
usage: patchOligos.py [-o org] [-c col] [-s cfg] infile outfile
  org    = organism (currently one of hs, mm, mm9)
  column = column in infile containing gene IDs
  cfg    = index of configuration to use

Existing configurations for {} (last one is used by default):
""".format(org)
        DB.listSettings()
