CC = g++
LDFLAGS = -DVERSION=$(shell git describe --tags) -pthread
CFLAGS = -Wall -Wextra -Wfatal-errors -pedantic -std=c++11
OPT = -O3
# OPT = -ggdb

SOURCE_DIR = src
BUILD_DIR = build
BIN_DIR = ~/bin/
DEPEND_FILE = $(BUILD_DIR)/.depend

EXEC = MakeShape Shape2Zernike Zernike2Shape
MODS = iotools mesh moments triangle vec zernike
OTHER = Makefile

TARGETS = $(addprefix $(BUILD_DIR)/, $(EXEC))
OBJ = $(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(MODS)))
OBJ_DEP = $(OTHER) | $(BUILD_DIR)

.SUFFIXES:
.PHONY: all bin doc clean

all: $(TARGETS)

bin: all
	cp -f $(TARGETS) $(BIN_DIR)

doc:
	cd $(SOURCE_DIR) && doxygen
	cd $(BUILD_DIR)/latex/ && make
	ln -f -s $(BUILD_DIR)/latex/refman.pdf
	ln -f -s $(BUILD_DIR)/html/index.html

clean:
	rm -rf $(BUILD_DIR) refman.pdf index.html

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

$(TARGETS): $(BUILD_DIR)/%: $(SOURCE_DIR)/%.cpp $(OBJ)
	$(CC) $(CFLAGS) $(OPT) -o $@ $< $(OBJ) $(LDFLAGS)

$(OBJ): $(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(OBJ_DEP)
	$(CC) $(CFLAGS) $(OPT) -o $@ -c $<

$(DEPEND_FILE): | $(BUILD_DIR)
	$(info creating $(DEPEND_FILE))
	$(foreach mod, $(MODS), $(eval($(shell g++ -MM -MT $(BUILD_DIR)/$(mod).o $(SOURCE_DIR)/$(mod).cpp >> $(DEPEND_FILE)))))
	$(foreach prog, $(EXEC), $(eval($(shell g++ -MM -MT $(BUILD_DIR)/$(prog) $(SOURCE_DIR)/$(prog).cpp >> $(DEPEND_FILE)))))

include $(DEPEND_FILE)
