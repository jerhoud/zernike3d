CC = g++
LDFLAGS = -DVERSION=$(shell git describe)
CFLAGS = -Wall -Wextra -Wfatal-errors -pedantic -std=c++11 -DZER_N200
OPT = -O3
# OPT = -ggdb

SOURCE_DIR = src
BUILD_DIR = build
BIN_DIR = ~/bin/

WITNESS = $(BUILD_DIR)/witness
EXEC = makeOFF zm rzm
MODS = iotools mesh moments triangle vec zernike_data zernike
HEADERS = zernike_def arg_parse
OTHER = Makefile

TARGETS = $(addprefix $(BUILD_DIR)/, $(EXEC))
OBJ = $(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(MODS)))
DEP = $(addsuffix .hpp, $(addprefix $(SOURCE_DIR)/, $(MODS) $(HEADERS))) $(OTHER) $(WITNESS)

all: $(TARGETS)


$(WITNESS):
	mkdir -p $(BUILD_DIR)
	touch $(WITNESS)

bin: all
	cp -f $(TARGETS) $(BIN_DIR)

$(TARGETS): $(BUILD_DIR)/%: $(SOURCE_DIR)/%.cpp $(OBJ) $(DEP)
	$(CC) $(CFLAGS) $(OPT) -o $@ $< $(OBJ) $(LDFLAGS)

doc:
	cd $(SOURCE_DIR) && doxygen
	cd $(BUILD_DIR)/latex/ && make
	ln -f -s $(BUILD_DIR)/latex/refman.pdf
	ln -f -s $(BUILD_DIR)/html/index.html

$(OBJ): $(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(DEP)
	$(CC) $(CFLAGS) $(OPT) -o $@ -c $<

clean:
	rm -rf build refman.pdf index.html
