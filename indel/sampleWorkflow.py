vtbl_Dict = {}
genos_Dict = {}
phased_vtbl_Dict = {}
phased_genos_Dict = {}

autosomalList = ["2L","2R","3L","3R"]

for name in autosomalList:
    
    chromChunk = indel.Chrom(name,'../data/variants/cross_indels.Aug29.chromgroup.h5')
    chromChunk.read_metadata("/Users/beccalove 1/Desktop/Indels/AG1K_indels/metadata/cross.txt")
    chromChunk.ID_positions(is_CDS)
    chromChunk.extract_exonic()
    chromChunk.filter_GQ_MQ_QD()
    chromChunk.filter_on_parents()
    chromChunk.ID_autosomal_Mendelian_violations()
    chromChunk.remove_Mendelian_violations()
    chromChunk.phase_and_filter(permitted_nonphased=len(chromChunk.metadata))

    vtbl_Dict[name] = chromChunk.vt_Mendel_filtered
    genos_Dict[name] = chromChunk.gt_Mendel_filtered
    phased_vtbl_Dict[name] = chromChunk.vt_phased
    phased_genos_Dict[name] = chromChunk.phased_genos


