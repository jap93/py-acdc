#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>


#include <cstdint>
#include <filesystem>

#ifdef __cpp_lib_filesystem
namespace fs = std::filesystem;
#endif
#ifdef __cpp_lib_experimental_filesystem
namespace fs = std::experimental::filesystem;
#endif

/**
 * @brief splits a string into a vector so that it is similar to python
 * 
 * @param s (std::string) : the string to be split
 * @return std::vector<std::string> : vector of component strings
 */
std::vector<std::string> split(std::string s)
{
    std::vector<std::string> words;

    std::stringstream ss(s);
    std::string word;
    
    while (ss >> word) 
    {
        words.push_back(word);
    }

    return words;
}

//adapted from the c++ reference
bool check_exists(const fs::path& p, fs::file_status s = fs::file_status{})
{
    //std::cout << p;
    if(fs::status_known(s) ? fs::exists(s) : fs::exists(p))
        return true;
    else
        return false;
}

void copyFile(std::string fileName, std::string destination)
{
    std::string destFile = destination + "/" + fileName;
    //first check the file exists
    fs::directory_entry tmp_file{fileName};
 
    if (tmp_file.exists()) 
    {
        std::string cmd = "cp " + fileName + " " + destFile;
        //copy to deistantion folder
        system(cmd.c_str());
    }
}

void dcell(const double* aaa, double* bbb)

/*
!
! dl_poly_3 subroutine to calculate the dimensional properies of a
! simulation cell specified by the input 3x3 matrix aaa [cell vectors in
! rows, the matrix is in the form of one dimensional reading
! [row1,row2,row3].
!
! The results are returned in the array bbb, with:
!
! bbb[1 to 3] - lengths of cell vectors
! bbb[4 to 6] - cosines of cell angles
! bbb[7 to 9] - perpendicular cell widths
! bbb[10]     - cell volume
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
{
  double
      axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3; 
      //x[3],y[3],z[3],d[3];

//    calculate lengths of cell vectors

      bbb[0] = sqrt(aaa[0]*aaa[0] + aaa[1]*aaa[1] + aaa[2]*aaa[2]);
      bbb[1] = sqrt(aaa[3]*aaa[3] + aaa[4]*aaa[4] + aaa[5]*aaa[5]);
      bbb[2] = sqrt(aaa[6]*aaa[6] + aaa[7]*aaa[7] + aaa[8]*aaa[8]);

//    calculate cosines of cell angles

      bbb[3]=(aaa[0]*aaa[3]+aaa[1]*aaa[4]+aaa[2]*aaa[5])/(bbb[0]*bbb[1]);
      bbb[4]=(aaa[0]*aaa[6]+aaa[1]*aaa[7]+aaa[2]*aaa[8])/(bbb[0]*bbb[2]);
      bbb[5]=(aaa[3]*aaa[6]+aaa[4]*aaa[7]+aaa[5]*aaa[8])/(bbb[1]*bbb[2]);

      //bbb[4] = acos(bbb[4]) / TORADIANS
      //bbb[5] = acos(bbb[5]) / TORADIANS
      //bbb[6] = acos(bbb[6]) / TORADIANS

//    calculate vector products of cell vectors

      axb1=aaa[1]*aaa[5]-aaa[2]*aaa[4];
      axb2=aaa[2]*aaa[3]-aaa[0]*aaa[5];
      axb3=aaa[0]*aaa[4]-aaa[1]*aaa[3];   

      bxc1=aaa[4]*aaa[8]-aaa[5]*aaa[7];
      bxc2=aaa[5]*aaa[6]-aaa[3]*aaa[8];
      bxc3=aaa[3]*aaa[7]-aaa[4]*aaa[6];
      
      cxa1=aaa[7]*aaa[2]-aaa[1]*aaa[8];
      cxa2=aaa[0]*aaa[8]-aaa[2]*aaa[6];
      cxa3=aaa[1]*aaa[6]-aaa[0]*aaa[7];
      
//    calculate volume of cell

      bbb[9]=abs(aaa[0]*bxc1+aaa[1]*bxc2+aaa[2]*bxc3);

//    calculate cell perpendicular widths

      bbb[6]=bbb[9]/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3);
      bbb[7]=bbb[9]/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3);
      bbb[8]=bbb[9]/sqrt(axb1*axb1+axb2*axb2+axb3*axb3);

}
