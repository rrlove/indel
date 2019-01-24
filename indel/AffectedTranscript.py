class AffectedTranscript():
    
    def __init__(self, chrom, name):
        
        self.chrom = chrom
        self.genos = None
        self.indices = []
        self.name = name
        self.positions = []
        self.ranges = []
        self.vtbl = None
        
    @property
    def n_exons(self):
        
        return len(self.ranges)
        
    def extract_vtbl(self, vtbl):
        
        if not len(self.indices) > 0:
            
            raise ValueError("This transcript has no associated variants")
            
        else:
            
            self.indices.sort()
            
            self.vtbl = vtbl[vtbl["CHROM"] ==\
                             self.chrom.encode('ascii')][[self.indices]]
            
    def extract_genos(self, genos):
        
        if not len(self.indices) > 0:
            
            raise ValueError("This transcript has no associated variants")
            
        else:
            
            self.indices.sort()
            
            self.genos = genos[[self.indices]]
