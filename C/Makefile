FLAGS=-O3 -std=gnu99
EXEC=knn
CC=gcc

ALL: $(EXEC)

$(EXEC):
	$(CC) $(FLAGS) $(EXEC).c -c -o $(EXEC).o
	$(CC) $(FLAGS) $(EXEC).o  -o $(EXEC)

clean:
	rm -rf $(EXEC) *.o
