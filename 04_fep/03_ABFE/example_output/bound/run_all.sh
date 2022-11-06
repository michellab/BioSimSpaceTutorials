for i in {1..5};
do
	cd run00$i;
	./runme.sh;
	cd ..;
done
