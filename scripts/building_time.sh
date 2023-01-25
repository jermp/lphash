/usr/bin/time ./lphash build-p /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz 63 18 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time ./lphash build-p /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz 63 20 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time ./lphash build-p /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz 63 24 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time ./lphash build-p /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz 63 24 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time ./lphash build-p /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz 63 28 -t 4 -c 5.0 2>> lphash_building_time.log

/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log

/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz 63 -b out.bin -t 4 -d /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log

/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz 63 -p out.bin -t 4 -c 5.0 -d /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz 63 -p out.bin -t 4 -c 5.0 -d /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz 63 -p out.bin -t 4 -c 5.0 -d /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz 63 -p out.bin -t 4 -c 5.0 -d /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time ./ptbb_build /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz 63 -p out.bin -t 4 -c 5.0 -d /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
