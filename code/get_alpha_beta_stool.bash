SHARED=$1 #data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.shared
SUBSAMPLE=1000

mothur "#summary.single(shared=$SHARED, calc=nseqs-coverage-invsimpson-shannon-sobs, subsample=$SUBSAMPLE);
				dist.shared(calc=braycurtis, subsample=$SUBSAMPLE)"

rm $(echo $SHARED | sed -e "s/shared/groups.summary/")
rm $(echo $SHARED | sed -e "s/shared/braycurtis.0.03.lt.dist/")
rm $(echo $SHARED | sed -e "s/shared/braycurtis.0.03.lt.std.dist/")
