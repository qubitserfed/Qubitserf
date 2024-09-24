import os
import sys

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "dev":
        cflags = "-std=c++17 -g -pthread"
    else:
        cflags = "-std=c++17 -Ofast -pthread"

    # build object files
    source_files = []
    library_objs = []
    for file in os.listdir('src'):
        file_name, file_ext = os.path.splitext(file)

        if os.path.isfile(os.path.join('src', file)) and file_ext=='.cpp':
            os.system(f'g++ -c {cflags} src/{file} -o build/{file_name}.o')
            if not file_name in ['testing', 'interface', 'library']:
                library_objs.append(file_name + '.o')

    library_deps = " ".join(map(lambda filename: "build/" + filename, library_objs))
    # build interface
    os.system(f'g++ {cflags} {library_deps} build/interface.o -o build/interface')

    # build tester
    os.system(f'g++ {cflags} {library_deps} build/testing.o -o build/testing')

    # build c library
    os.system(f'g++ -c {cflags} {library_deps} build/library.o -o build/library')
