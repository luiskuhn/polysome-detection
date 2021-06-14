from polysomeDetection.basicStructures import PolysomeParticleList2, PolysomeParticleList



pl_path = '../../../datasets/inVitro/poly_cyto_1_aligned_full_bin1_mirr.xml'

print pl_path

pl = PolysomeParticleList()
pl.fromXMLFile(pl_path)
print 'particles', len(pl)

pl2 = PolysomeParticleList2(pl=pl, polysomeType='testType')
print 'particles', len(pl2)

pl2.toXMLFile('test.xml')

