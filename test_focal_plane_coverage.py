# Series of tests for focal_plane coverage
import parameters as par
import numpy as np

def test_pitch(pitch, test_pitch):
    tol = 1e-3
    if pitch - test_pitch <= tol:
        print('Pitch test passed')
        print(test_pitch)
    else:
        
        print('Pitch test failed: review x_inc and y_inc definition')

def main():
    test_pitch(par.pitch, par.test_pitch)

if __name__ == '__main__':
    main()