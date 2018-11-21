SHARED=$1 #data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.shared
SUBSAMPLE=1000

cut -f 2 $SHARED | grep "DNA" | sed -E "s/((.*)_DNA.)/\1\t\2/" > data/mothur/merge.file

mothur "#merge.groups(shared=$SHARED, design=data/mothur/merge.file);
				summary.single(calc=nseqs-coverage-invsimpson-shannon-sobs, subsample=$SUBSAMPLE);
				dist.shared(calc=braycurtis, subsample=$SUBSAMPLE);
				dist.shared(shared=$SHARED, calc=braycurtis, subsample=500)"

rm data/mothur/merge.file
rm $(echo $SHARED | sed -e "s/shared/merge.shared/")
rm $(echo $SHARED | sed -e "s/shared/merge.groups.summary/")
rm $(echo $SHARED | sed -e "s/shared/merge.braycurtis.0.03.lt.dist/")
rm $(echo $SHARED | sed -e "s/shared/merge.braycurtis.0.03.lt.std.dist/")
rm $(echo $SHARED | sed -e "s/shared/braycurtis.0.03.lt.dist/")
rm $(echo $SHARED | sed -e "s/shared/braycurtis.0.03.lt.std.dist/")
