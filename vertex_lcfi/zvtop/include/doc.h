
/*! \page ZVTOP


 The main algorithm for ZVRES is VertexFinderClassic  and for ZVKIN VertexFinderGhost
 
 "classic" denotes implementation as in Dave Jacksons original ZVTOP paper - Nuc. Inst. Meth. A 388 (1997) 247-253
 
 Please note: This documentation is not yet complete - it needs references to classes
 
 
 
 
ZVTOP is a proven vertex finding algorithm developed by D Jackson and used at SLD.  It consists of two complementary vertex finding methods, ZVTOP and ZVKIN which are detailed below. ZVTOP was chosen for this study due to its previous use for ILC studies and the local presence of the initial developer of the algorithm. The ability of the ZVKIN branch to reconstruct decays with only one seen track was also seen as desirable for the charge dipole technique.

\subsection{ZVRES}

ZVRES is the topological branch of ZVTOP, most suited to decays where there is more than one seen track from each vertex. The algorithm proceeds by locating points in space where vertices are likely to exist through use of a heuristic function that measures track density at a given point. Ambiguities in which tracks are assigned to which vertices are then resolved by comparing the magnitude of this function.

Where the version of the algorithm presented here differs from the original ZVTOP paper this is due to the paper differing from the SGV FORTRAN implementation that has been used for previous studies. To maintain performance and validity of comparisons the FORTRAN implementation is considered authoritative.

Before detailing the steps of the algorithm it is necessary to define some devices used. The effect of parameters introduced is discussed after the detail of the algorithm.

\textit{\textbf{Vertex Function}}\\
 \f$ V(\mathbf{r}) \f$  is the heuristic function used to measure track density at a given point. It is based on a Gaussian tube function  \f$ f_i(\mathbf{r}) \f$ ; in the original ZVTOP paper it is defined for each track:
\f[
f_i(\mathbf{r})=\exp\lbrace-\frac{1}{2}[(\frac{x'-(x'_0+\kappa y'^2)}{\sigma_T})^2+(\frac{z-(z_0+\tan(\lambda)y')}{\sigma_L})^2]\rbrace
\f]  

Note that this is a parabolic approximation to the track trajectory where  \f$ x' \f$  and  \f$ y' \f$  are such that the track momentum is parallel to the  \f$ y' \f$  axis at the track's point of closest approach to the IP.  \f$ x'_0 \f$  and  \f$ z_0 \f$  are the co-ordinates at this point.  \f$ \kappa \f$  is the curvature in the  \f$ xy \f$  plane and  \f$ \lambda \f$  the track's dip angle in the  \f$ yz \f$  plane with respect to the y axis.  \f$ \sigma_L \f$  and  \f$ \sigma_{T} \f$  are the track errors in the  \f$ z \f$  direction and  \f$ xy \f$  plane respectively. For the version of ZVTOP developed for these studies the parabolic approximation was removed and replaced with:
\f[
f_i(\mathbf{r})=\exp\lbrace-\frac{1}{2}(\mathbf{r}-\mathbf{p})\mathbf{V}^{-1}(\mathbf{r}-\mathbf{p})^T\rbrace
\f]  

Where vector  \f$ \mathbf{p} \f$  is the point of closest approach on the track to  \f$ \mathbf{r} \f$  and matrix  \f$ \mathbf{V} \f$  is the 3D covariance of the track at  \f$ \mathbf{p} \f$ .  \f$ f_i(\mathbf{r}) \f$  is a spacial function that has a Gaussian cross section proportional to the covariance of the track at any given point. Note that the Gaussian is unnormalised so that the function is the weighted distance to the track relative to all other tracks.

These independent track functions are combined to produce a vertex function that gives a measure of the likelyhood of a track co-incidence at any point:
\f[
V(\mathbf{r})=\sum^N_{i=1}f_i(\mathbf{r})-\frac{\sum^N_{i=1}f^2_i(\mathbf{r})}{\sum^N_{i=1}f_i(\mathbf{r})} 
\f]  

Note that this function is near-zero at points that are near zero or one tracks, and has peaks near the coincidences of two or more tracks and can therefore be used to identify points in detector space that contain vertices. This represents the most basic form of  \f$ V(\mathbf{r}) \f$ , it is modified to suppress fake vertices by two mechanisms: the addition of an IP object and a weighting by distance from the jet axis.

The IP object is added by defining  \f$ f_0(\mathbf{r}) \f$  for the IP analogously to  \f$ f_i(\mathbf{r}) \f$  for the tracks:
\f[
f_0(\mathbf{r})=\exp\lbrace-\frac{1}{2}(\mathbf{r}-\mathbf{p})\mathbf{V}^{-1}(\mathbf{r}-\mathbf{p})^T\rbrace
\f]  



Where  \f$ \mathbf{p} \f$  is the IP position and  \f$ \mathbf{V} \f$  the IP covariance. This  \f$ f_0(\mathbf{r}) \f$  is inserted into  \f$ V(\mathbf{r}) \f$  with a weight  \f$ K_{IP} \f$ :
\f[
V(\mathbf{r})=K_{IP}f_0(\mathbf{r})+\sum^N_{i=1}f_i(\mathbf{r})-\frac{K_{IP}f^2_0(\mathbf{r})+\sum^N_{i=1}f^2_i(\mathbf{r})}{K_{IP}f_0(\mathbf{r})\sum^N_{i=1}f_i(\mathbf{r})} 
\f]  

This ensures a large peak in  \f$ V(\mathbf{r}) \f$  at the IP location, which serves to suppress fake vertices near the IP.

A weighting is implemented to suppress vertices that are far removed from the jet axis which are likely to be fake. This takes the form of a cylinder around the jet axis in which  \f$ V(\mathbf{r}) \f$  is unchanged and outside of which  \f$ V(\mathbf{r}) \f$  decays in proportion to the angle  \f$ \alpha \f$  between the side of the cylinder and the line from the base of the cylinder to  \f$ \mathbf{r} \f$  (FIGURE):
\f[
V(\mathbf{r})\rightarrow V(\mathbf{r})\exp(-K_\alpha\alpha^2)
\f]  

The parameter  \f$ K_\alpha \f$  allows the divergence of the cone in which V(r) is reduced to be controlled

The final  \f$ V(\mathbf{r}) \f$  used is the one including IP and jet axis mechanisms, but the algorithm can be made to work without these. Any given vertex has two values of  \f$ V(\mathbf{r}) \f$ : 
 \f$ V(\mathbf{r_{VERT}}) \f$  measured at the fitted vertex position  \f$ V(\mathbf{r_{MAX}}) \f$  the local maximum in  \f$ V(\mathbf{r}) \f$  found by gradient descent from the fitted vertex position.

\textit{\textbf{Vertex Resolution}}\\
Another important definition is that of resolvability. This is performed using  \f$ V(\mathbf{r}) \f$ . Two vertices are resolved if:
\f[
\frac{\min\lbrace V(\mathbf{r}): \mathbf{r}\in \mathbf{r_1} + \alpha(\mathbf{r_2} - \mathbf{r_1}), 0\leq\alpha\leq1\rbrace}{\min\lbrace V(\mathbf{r_1}),V(\mathbf{r_2})\rbrace}<R_0
\f]  

Where  \f$ R_0 \f$  is a threshold parameter. In other words vertices are resolved if the ratio of the minimum  \f$ V(\mathbf{r}) \f$  between the vertices to the lowest  \f$ V(\mathbf{r}) \f$  at either vertex is lower than the threshold. In this way if the vertices have no significant valley between them they are not resolved. FIGURE

When two vertices are to be merged the one with the smallest  \f$ V(\mathbf{r_{MAX}}) \f$  is destroyed and its tracks and IP (if any) added to the other.

\textit{\textbf{Algorithm}}\\
The algorithm proceeds as follows:
Initially all possible track pairs are formed from tracks that are contained in the input jet. Each of these pairs are fitted to a common vertex. The pairs are then cut if:

Either tracks contribution to the fitted vertex  \f$ \chi^2 \f$  is greater than parameter  \f$ \chi^2_0 \f$ 

The value of  \f$ V(\mathbf{r}) \f$  at the fitted vertex position  \f$ 
(V(\mathbf{r_{VERT}})) \f$  is less than  \f$ V_0 \f$ 

If information about the IP is to be used a further set of pairs is created consisting of the IP and one jet track each. Again these are fitted and cut away if and only if the track or IP's contribution to the fitted vertex  \f$ \chi^2 \f$  is higher than  \f$ \chi^2_0 \f$ .

Note that for a two object fit such as this the  \f$ \chi^2 \f$  of both objects in the fit is equal so this crierion is equal to a cut on vertex  \f$ \chi^2 \f$  of  \f$ 2\chi^2_0 \f$ .

Each of the track-track or track-IP pairs remaining is a candidate vertex. A track is therefore associated with many candidate vertices at this stage. Each track is removed from all candidate vertices in an unresolved set apart from the one with highest  \f$ V(\mathbf{r_{VERT}}) \f$  in the set using the following method:

First for each track a set  \f$ CV_i \f$  is created of all candidate vertices that include track  \f$ i \f$ . For each vertex  \f$ X \f$  of  \f$ CV_i \f$  that has a  \f$ V(\mathbf{r_{VERT}}) \f$  less than  \f$ 10\% \f$  of the maximum  \f$ V(\mathbf{r_{VERT}}) \f$  in  \f$ CV_i \f$  track  \f$ i \f$  is removed from  \f$ X \f$  and that vertex removed from  \f$ CV_i \f$ .

The next step is repeated until  \f$ CV_i \f$  is empty: The vertex in  \f$ CV_i \f$  with the highest  \f$ V(\mathbf{r_{VERT}}) \f$  is removed from  \f$ CV_i \f$  and added to set  \f$ RV_i \f$ . Every vertex of  \f$ CV_i \f$  that is unresolved from any vertex in  \f$ RV_i \f$  has track  \f$ i \f$  removed and is removed from  \f$ CV_i \f$ . Note that the vertices are \textit{not} refit when these tracks are removed. They retain their inital two-track fitted position.

The IP is removed from all vertices except that with the highest  \f$ V(\mathbf{r_{VERT}}) \f$  that contains the IP. At this point all candidate vertices that have 0 objects (track or ip) are removed. Note that this leaves candidate vertices that have only one track, but which still retain their position from the inital two-track fit.

Each unresolved set of candidate vertices is now merged into a single vertex. This is performed by taking a candidate vertex and finding all vertices unresolved from it, then finding all vertices unresolved from them until no more are found to be unresolved. All vertices found in this unresolved set are then merged. (Diagram in terms of digraph)

Tracks are then cut from the vertices based on their  \f$ \chi^2 \f$  contribution.
If the track with the highest  \f$ \chi^2 \f$  contribution has a  \f$ \chi^2 \f$  contribution above the cut threshold, it is removed and the vertex refit. This is repeated till the track with the highest  \f$ \chi^2 \f$  contribution is below the cut threshold ( \f$ \chi^2_{0TRIM} \f$ )or the vertex has less than two tracks (less than one if the vertex has the ip object). After this  \f$ \chi^2 \f$  trimming all candidate vertices that are no longer define points in space are removed - ie. those with two or more tracks or that have the ip object are retained.

The last stage is to remove any remaining ambiguities in track allocation where a track is contained in more than one vertex. Ambiguities are resolved by retaining the track in the vertex with highest  \f$ V(\mathbf{r_{MAX}}) \f$  out of all vertices that contain that track. In descending order of  \f$ V(\mathbf{r_{MAX}}) \f$  the candidate vertices tracks are removed from all candidate vertices with a smaller value of  \f$ V(\mathbf{r_{MAX}}) \f$ .

\textit{\textbf{Parameters}}\\
 \f$ K_{IP} \f$  is the weight of the IP object in the vertex function. Due to the vertex function being used to resolve vertices an increased value of  \f$ K_{IP} \f$  will merge vertices near the IP with the IP. This can be useful for suppression of fake vertices, but note that this should not be necessary if the IP's covariance is correct as a larger covariance has the same effect as a higher  \f$ K_{IP} \f$ . The default value is 1.

 \f$ K_\alpha \f$  controls the divergence of the cone in which  \f$ V(\mathbf{r}) \f$  is non-zero. A larger value makes a less divergent cone. As lower energy jets are less collimated this value is made proportional to the jet energy. The default value is 0.125 times the jet energy in GeV.

 \f$ R_0 \f$  is the resolution criterion. A higher value will lead to more merging of vertices that are spatially close The default value is 0.6.
 \f$ \chi^2_0 \f$  and  \f$ V_0 \f$  determine which of the  \f$ \frac{1}{2}N(N+1) \f$  possible two object fits are immediately discarded at the beginning of the algorithm a higher value of  \f$ \chi^2_0 \f$  or lower  \f$ V_0 \f$  increases the number considered but will lead to fake vertices. The default values are 10 and 0.001 respectively.

 \f$ \chi^2_{0TRIM} \f$  controls which tracks are rejected before the final ambiguity resolution, a higher value reduces the number of tracks rejected. The default value is 10.

*/


