# WassersteinFisherRaoDistance
An optimization method for computing the unbalanced Wasserstein-Fisher-Rao optimal transport distance between two measures on S^2. Includes an application to computing the SRNF shape distance and color tranfer.

# References

```
@article{bauer2022SRNF, 
      author = {Martin Bauer and Emmanuel Hartman and Eric Klassen}, 
      title = {The Square Root Normal Field Distance and Unbalanced Optimal Transport}, 
      journal = {Applied Mathematics {\&} Optimization},  
      doi = {10.1007/s00245-022-09867-y}, 
      url = {https://link.springer.com/article/10.1007/s00245-022-09867-y}, 
      year = 2022, month = {may}, publisher = {Springer Science and Business Media {LLC}}, volume = {85}, number = {3}}
```
Please cite this paper in your work. Available at https://arxiv.org/abs/2105.06510

# Requirements

Our code is written in Python and includes the following Python libraries:

    -Numpy 1.19.5

    -Scipy 1.4.1

    -PyTorch 1.4.0

    -MatPlotLib 3.2.0
  
  For application to the SRNF shape distance:
  
    -open3d  0.10.0.1
  
    -polyhedrec (Available at https://github.com/gsellaroli/polyhedrec)
  
  For application to Color Transfer:
  
    -scikit-image 0.15.0

# Usage

Usage of WFR explained a jupyter notebook called *Examples.ipynb*. 

Usage of WFR for SRNF shape distance explained in jupyter notebook called *ShapeExample.ipynb*. 

Usage of WFR for color transfer explained in jupyter notebook called *ColorTransferExample.ipynb*.

# Main Functions

Discriptions of the main functions can be found by calling help(WFR),help(PLShapes), or help(ColorTransfer).
        
# License

You may redistribute and/or modify this code under the terms of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

# Contact

elh18e(at)my.fsu.edu

bauer(at)math.fsu.edu
