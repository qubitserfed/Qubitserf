# Compiler
CXX = g++

# Directories
SRC_DIR = src
BUILD_DIR = build

# Source files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)

# Object files (replace src/ with build/ and .cpp with .o)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))

# Flags
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -ffast-math -pthread -fPIC

# Targets
all: $(OBJ_FILES)

# Rule to build object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

OBJ_FILES_NO_INTERFACE_NO_TEST := $(filter-out $(BUILD_DIR)/interface.o $(BUILD_DIR)/testing.o, $(OBJ_FILES))

# Compile interface.cpp into an executable, linking with all object files except testing.o
interface: $(BUILD_DIR)/interface.o $(OBJ_FILES_NO_INTERFACE)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/interface $(BUILD_DIR)/interface.o $(OBJ_FILES_NO_INTERFACE_NO_TEST)

# Compile testing.cpp into an executable, linking with all object files except interface.o
testing: $(BUILD_DIR)/testing.o $(OBJ_FILES_NO_testing)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/testing $(BUILD_DIR)/testing.o $(OBJ_FILES_NO_INTERFACE_NO_TEST)

# Clean build directory
clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/interface $(BUILD_DIR)/testing
