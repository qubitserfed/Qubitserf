# Compiler
CXX = g++

# Directories
SRC_DIR = src
BUILD_DIR = build
TEST_DIR = tests

# Source files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
TEST_SRC_FILES := $(wildcard $(TEST_DIR)/test_*.cpp)
TEST_MAIN_SRC := $(TEST_DIR)/main.cpp

# Object files (replace src/ with build/ and .cpp with .o)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))
TEST_OBJ_FILES := $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/tests/%.o, $(TEST_SRC_FILES))
TEST_MAIN_OBJ := $(BUILD_DIR)/tests/main.o

# Flags
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -ffast-math -pthread -fPIC

# Targets
all: $(OBJ_FILES)

# Rule to build object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to build test object files
$(BUILD_DIR)/tests/%.o: $(TEST_DIR)/%.cpp | $(BUILD_DIR)/tests
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/tests:
	mkdir -p $(BUILD_DIR)/tests

OBJ_FILES_NO_INTERFACE_NO_TEST := $(filter-out $(BUILD_DIR)/interface.o $(BUILD_DIR)/testing.o, $(OBJ_FILES))

# Compile interface.cpp into an executable, linking with all object files except testing.o
interface: $(BUILD_DIR)/interface.o $(OBJ_FILES_NO_INTERFACE_NO_TEST)
	$(CXX) -Wall -Wextra -std=c++17 -O3 -ffast-math -pthread -o $(BUILD_DIR)/interface $(BUILD_DIR)/interface.o $(OBJ_FILES_NO_INTERFACE_NO_TEST)

# Compile testing.cpp into an executable, linking with all object files except interface.o
testing: $(BUILD_DIR)/testing.o $(OBJ_FILES_NO_INTERFACE_NO_TEST)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/testing $(BUILD_DIR)/testing.o $(OBJ_FILES_NO_INTERFACE_NO_TEST)

unit_tests: $(TEST_MAIN_OBJ) $(TEST_OBJ_FILES) $(OBJ_FILES_NO_INTERFACE_NO_TEST)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/unit_tests $(TEST_MAIN_OBJ) $(TEST_OBJ_FILES) $(OBJ_FILES_NO_INTERFACE_NO_TEST)

test: unit_tests
	./$(BUILD_DIR)/unit_tests

# Clean build directory
clean:
	rm -rf $(BUILD_DIR)

# Debug build
debug:
	$(MAKE) unit_tests BUILD_DIR=build/debug CXXFLAGS="-g -O0 -Wall -Wextra -std=c++17 -pthread -fPIC"
 