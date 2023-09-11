THIS_FULL_PATH=$(dirname $(debian_realpath "$0"))
LPHASH_DIR=$(dirname "$THIS_FULL_PATH")
LPBUILD="$LPHASH_DIR/build/lphash partitioned build"
PTBUILD="$LPHASH_DIR/build/ptbb_build"

/usr/bin/time $LPBUILD -i /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz -k 63 -m 18 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time $LPBUILD -i /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz -k 63 -m 20 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time $LPBUILD -i /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time $LPBUILD -i /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -m 24 -t 4 -c 3.0 2>> lphash_building_time.log
/usr/bin/time $LPBUILD -i /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -m 28 -t 4 -c 5.0 2>> lphash_building_time.log

/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 2 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -b out.bin -t 4 -d -i /data2/DNA/tmp_dir/ -g 1 2>> bbhash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/yeast.k63.unitigs.fa.ust.fa.gz -k 63 -p out.bin -t 4 -c 5.0 -d -i /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/celegans.k63.unitigs.fa.ust.fa.gz -k 63 -p out.bin -t 4 -c 5.0 -d -i /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/cod.k63.unitigs.fa.ust.fa.gz -k 63 -p out.bin -t 4 -c 5.0 -d -i /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/kestrel.k63.unitigs.fa.ust.fa.gz -k 63 -p out.bin -t 4 -c 5.0 -d -i /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
/usr/bin/time $PTBUILD -i /data2/DNA/lphash_datasets/human.k63.unitigs.fa.ust.fa.gz -k 63 -p out.bin -t 4 -c 5.0 -d -i /data2/DNA/tmp_dir/ 2>> pthash-v2_building_time.log
