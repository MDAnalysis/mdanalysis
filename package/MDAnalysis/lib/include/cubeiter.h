#ifndef DISTOPIA_CUBEITER_H
#define DISTOPIA_CUBEITER_H

/*
 * usage:
 * ConcentricCubeIterator ci();
 *
 * ci.x, ci.y, ci.z  // all 0,0,0 at this point
 * ci.next();  // will return false, only one cube in central shell
 * ci.expand();  // moves onto 3x3 shell of cubes
 * while (ci) {
 *   ci.x, ci.y, ci.z  // coordinate of a cube in concentric shell
 *   ci.next();  // moves onto next cube in this concentric shell
 * }
 * ci.expand()  // moves onto 5x5 shell of cubes etc
 */
struct ConcentricCubeIterator {
    int x, y, z;
    unsigned short shell_number;
    char face;  // [0-5] for which face we're currently on, -1 for no face

    ConcentricCubeIterator() : x(0), y(0), z(0), face(0), shell_number(0) {};

    // advance iterator to next position, returns if now valid to access x, y, z
    bool next() {
      if (shell_number == 0) {
        face = -1;
        return false;
      }
      switch (face) {
        /*
         * 0: -x, yz
         * 1: +x, yz
         * 2: -y, xz (limit x bounds by 1 each side)
         * 3: +y, xz
         * 4: -z, xy (limit x & y bounds by 1 each side)
         * 5: +z, xy
         */
        case 0:
        case 1:
          y++;
          if (y == shell_number + 1) {
            y = - shell_number;
            z++;
            if (z == shell_number + 1) {
              if (face == 0) {
                x = shell_number;
                y = - shell_number;
                z = - shell_number;
              } else {
                x = - (shell_number - 1);
                y = - shell_number;
                z = - shell_number;
              }
              face++;
            }
          }
          break;
        case 2:
        case 3:
          x++;
          if (x == shell_number) {
            x = - (shell_number - 1);
            z++;
            if (z == shell_number + 1) {
              if (face == 2) {
                x = - (shell_number - 1);
                y = shell_number;
                z = - shell_number;
              } else {
                x = - (shell_number - 1);
                y = - (shell_number - 1);
                z = - shell_number;
              }
              face++;
            }
          }
          break;
        case 4:
        case 5:
          x++;
          if (x == shell_number) {
            x = - (shell_number - 1);
            y++;
            if (y == shell_number) {
              if (face == 4) {
                x = - (shell_number - 1);
                y = - (shell_number - 1);
                z = shell_number;
                face++;
              } else {
                x = 0;
                y = 0;
                z = 0;
                face = -1;
              }
            }
          }
          break;
        default:
          break;
      }
      return face != -1;
    }

    operator bool() {
      return face != -1;
    }

    // move iterator onto next concentric shell
    void expand() {
      shell_number++;

      face = 0;
      x = - shell_number;
      y = - shell_number;
      z = - shell_number;
    };

    void reset() {
      shell_number = 0;
      face = 0;
      x = 0;
      y = 0;
      z = 0;
    }
};

#endif //DISTOPIA_CUBEITER_H
