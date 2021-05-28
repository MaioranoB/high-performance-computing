import subprocess

def make_files():
    # -fopenmp flag in all calls to use omp_get_wtime()
    subprocess.run('g++ -fopenmp -o laplace.out laplace.cxx'.split())
    subprocess.run('g++ -fopenmp -O3 -o laplace_O3.out laplace.cxx'.split())
    subprocess.run('g++ -fopenmp -o laplace_fopenmp.out laplaceOpenMP.cxx'.split())
    subprocess.run('g++ -fopenmp -O3 -o laplace_fopenmp_O3.out laplaceOpenMP.cxx'.split())

def main():
    NX = [512,1024,2048]
    flags = ['', '_O3', '_fopenmp', '_fopenmp_O3']
    file = open('report.csv', 'w')

    for flag in flags:
        executable = ''.join(['./laplace', flag, '.out'])
        for nx in NX:
            command = [executable,str(nx)]
            output = subprocess.run(command, stdout = file)
            # handle(output.returncode)

    file.close()

if __name__ == '__main__':
    make_files()
    main()
