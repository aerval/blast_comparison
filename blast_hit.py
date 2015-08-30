from Bio import Entrez  # To access NCBIs Genebank database


class BlastHit(object):
    '''
    Object that contains all data about a particular Hit
    '''

    def __init__(self, hit):
        '''
        Notice: some of the properties included here do not appear in the
        short tabular BLAST representation. If one would adapt the program to
        also work with that kind of data, all the properties of higher hit
        indices (from hit[13]?) would be exluded. Note that you would also
        need to adjust compare_hit function to only check the primary
        properties (and further that this would mask hits with changes in the
        description as it is not included in the short tabular format).
        '''

        self.raw = hit
        hit = hit.split('\t')

        # Loop through the ID fields and extract the different IDs
        ids = hit[1].strip('|').split('|')
        self.ids = []
        for n, i in enumerate(ids):
            if len(i) < 5:
                # 5 is a relativly arbitrary number to distiguish between hit
                # ids and databank ids
                num = n + 1
                while num < len(ids):
                    if len(ids[num]) >= 5:  # Very vague, it might be better to
                        # have a list of potential dbs.
                        self.ids.append(GeneId(ids[n], ids[num]))
                        # -> Every id is a list containing the database and
                        # the id itself
                        num += 1
                    else:
                        break

        self.name = hit[0]
        self.length = hit[3]
        self.score = hit[13]
        self.e_value = hit[10]
        # self.bit_score = hit[11]
        # -> Problematic due to 12261 != 1.226e04.  This can appear probably
        # at many places but here i detected it taken out because very similar
        # to e-value
        self.identities = hit[14]
        self.positives = hit[15]
        self.gaps = hit[16]
        self.qframe = hit[18]
        self.sframe = hit[19]
        self.query = hit[20]
        self.query_start = hit[6]
        self.mismatch = hit[4]
        self.subject = hit[21]
        self.subject_start = hit[8]
        self.status = ''
        # States whether a similar BlastHit has been found in the other Search

    def compare_hit(self, other, check_ids=True):
        '''
        Returns true in those cases where only eValue (or description) may be
        deviating and every other stat is equal
        
        compare = the BlastHit the current hit is compared to
        check_ids = wether the comparison should also look for the same ids.
        It sometimes happens that a database entree is updated with new
        information while the basic properties do not changes (and therefore
        also not the hit charateristics). Check_ids = True therefore checks
        for Ids as well while = False ignores it.
        '''

        differences = []

        # if self.bit_score != other.bit_score: # Should be the same as score
        #    return False, None                 # -> see above
        if check_ids and not self.ids == other.ids:
            return False, None
        elif self.query_start != other.query_start:
            return False, None
        elif self.subject_start != other.subject_start:
            return False, None
        elif self.qframe != other.qframe:
            return False, None
        elif self.sframe != other.sframe:
            return False, None
        elif [letter for letter in self.query if letter in 'ACGT'] !=
             [letter for letter in other.query if letter in 'ACGT']:
            # Neccesary due to alignment differences like 'C-T' != 'CT-' which
            # come from changes in alignment algorithm.
            return False, None
        elif [letter for letter in self.subject if letter in 'ACGT'] !=
             [letter for letter in other.subject if letter in 'ACGT']:
            # see above "query"
            return False, None
        else:
            # Check for smaller differencens (these are database changes that
            # do not directly change the hit.  For example, when the database
            # becomes bigger, this make different (hit) cases more likely. The
            # reason why these are nevertheless checked is to clarify that
            # these are no perfect hits and therefore line by line comparison
            # would not find this (also, a slightly higher e-Value might lead
            # to exclusion from the BLAST algorithm).
            if self.e_value != other.e_value:
                differences.append('eValue')
                
            # Changes in alignment (due to internal changes in BLAST like
            # "CTGC" and "CATC" as gapless alignment or "C-TGC" and "CAT-C").
            if self.identities != other.identities or 
               self.query != other.query or
               self.subject != other.subject or
               self.score != other.score or
               self.positives != other.positives or
               self.gaps != other.gaps or
               self.mismatch != other.mismatch:
                differences.append('alignment')
        return True, differences
        
    def __str__(self):
        '''
        When a BLAST hit is called as a string return the orginal data and,
        if specified, the status.
        '''
        
        if self.status == '':
            return self.raw
        else:
            return '%s\t%s' % (self.raw[:-1], self.status)
    
class GeneId(object):
    '''
    Class that describes an ID as it is found in a certain database
    
    db = Identifier of the database (eg. 'gb')
    num = identification number of the gene in the database
    '''
    
    def __init__(self, db, num):
        self.db = db
        self.num = num

    def __eq__(self, id2):
        '''
        Comparison of two GeneIds
        '''
        
        return self.db == id2.db and self.num == id2.num
    
def compare_blasts(blastA, blastB):
    '''
    Compare two lists of BlastHits. Returns two dictionaries that devide each
    lists in the BlastHits that appear in both lists ('same'), that appear
    with slight changes (in description or eValue) in both lists ('similar')
    and the Hits that appear in only one list ('unknown')
    '''
    
    hitsA = {'same':[], 'similar':[], 'unknown':[], 'all':list(blastA)}
    # -> The all list is inserted as copy to prevent overwriting
    hitsB = {'same':[], 'similar':[], 'unknown':[], 'all':list(blastB)}
    
    for hitA in blastA:
        found = False
        for num, hitB in enumerate(blastB):
            same, differences = hitA.compare_hit(hitB)
            if same:
                found = True
                if len(differences) > 0:
                    hitA.status = 'similar'
                    hitB.status = 'similar'
                    hitsA['similar'].append(hitA)
                    hitsB['similar'].append(hitB)
                else:
                    hitA.status = 'equal'
                    hitB.status = 'equal'
                    hitsA['same'].append(hitA)
                    hitsB['same'].append(hitB)
                blastB.pop(num)
                break 
                # NOTICE: This assumes that there is only one matching hit in
                # the other BLAST Search.  As long as only (description and)
                # e-value count as possible differences this is acceptable
                # (because very unlikely since basically impossible) but might
                # get problematic as soon as other properties are seen as
                # optional
        
        if not found:
            hitsA['unknown'].append(hitA)
    hitsB['unknown'].extend(blastB)
    
    return hitsA, hitsB

def get_idlist(hits, email):
    '''
    For a given list of Genebank Ids, this function retruns (in form of a dict
    with ids as key) some basic data like creation and updating date,
    description and whether the entree is still valid.
    
    hits = list of BlastHits
    email = Email Address (string) for accessing Entrez. They requiere this to
    let you notice if they detect excess use of their service
    '''
    
    # isolate the genebank ids
    ids = []
    for hit in hits:
        for ID in hit.ids:
            if ID.db == 'gi':
                if not ID.num in ids:
                    ids.append(ID.num)
    
    # Access Entrez and download the basic info
    Entrez.email = email
    # Access the nucleotide database
    # If you have another type of data (like protein) change the database name
    handle = Entrez.esummary(db='nucleotide', id=','.join(ids))
    records = Entrez.read(handle)
    
    gene_summaries = {}
    
    # The records can be accessed in the gene_summaries dictionary via their
    # ids
    for record in records:
        gene_summaries[str(record['Gi'])] = record
        
    return gene_summaries
    