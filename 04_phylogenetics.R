###### Lesson on phylogenetics in R

###### Places to find help:

# Analyses of Phylogenetics and Evolution with R by Emmanuel Paradis (author of ape) can find ebook pdf on google

# R-sig-phylo email list serv: https://www.mail-archive.com/r-sig-phylo@r-project.org/b

# Phytools blog: Liam Revell (author of phytools): blog.phytools.org

###### Essential packages! Install if you haven't already

library(ape)   #the "original" one. Needed to read and plot trees in R
library(geiger)   #expands capabiltiies of ape, some macroevolutionary analyses
library(phytools)   #mostly for tree manipulation
library(caper)   #includes even more analyses like PGLS


##### Tree structure in R

text.string<-
  "(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark);"

vert.tree<-read.tree(text=text.string)

### Try printing the vert.tree object to the screen to get a sense for what the tree object looks like in R.

vert.tree
str(vert.tree)

## The tree object contains a set of lists: $edge, $tip.label, $Nnode

vert.tree$tip.label  # Is a list of tip labels in the tree.

## Now plot the tree. The simplest way is just using plot()

plot(vert.tree)

## You can make some changes to the appearance easily:

plot(vert.tree, edge.width = 2)
plot(vert.tree, edge.width = 2, type="fan")

### We'll explore the attributes of trees in R a little further. Tips on the tree are numbered from 1 through the total number. Then the nodes are numbered from there, starting at the root. This is useful knowledge for a lot of tree manipulation functions.

plot(vert.tree) 
nodelabels() #visualize the node nabels. So this tree has 11 tips. So the root node is #12, and nodes are counted out from there. The number of nodes will always be equal to #tips minus 1. 

### Let's work on basic tree manipulation. We want a tree of just the mammals. How do we do that? Two basic ways: cut out anything that isn't a mammal, or take away the clade of mammals.

new.tree <- drop.tip(vert.tree, c("shark", "gold_fish", "coelacanth", "iguana", "robin"))

## notice that we added the underscore for gold fish. If you go vert.tree$tip.label, you notice it is written as gold_fish in the list object.

plot(new.tree)

### Or, just extract the clade for mammals. This requires you to know what node number belongs to that clade. 

plot(vert.tree) 
nodelabels()  # So the mammal clade is node #16

new.tree2 <- extract.clade(vert.tree, 16)
plot(new.tree2) 

## new.tree and new.tree2 are identical.


##### Now let's work with a gigantic tree. This will force us to use some tree manipulation functions.

##### This is a tree of all ray-finned fishes from Rabosky et al. 2013 (over 7,000 species!!!!!!)

http://datadryad.org/bitstream/handle/10255/dryad.48798/Rabosky_et_al_timetree.tre?sequence=1

## When you open this page in your browser, you notice it follows the parenthetical stucture above like for the small tree, but there are numbers following the tip labels. Those are branch lengths. This is a time calibrated tree, so the branch lengths are in units of millions of years. 

## To download this tree, hit Select-All, then copy and paste into a text editor like Text Wrangler. Save in your working directory with the file extension ".tree"

### Do NOT plot this tree in the R console! It is big so it will take a lot of time to plot. We are going to work on cutting this tree down to something more useful.

# Read the tree from your file:

fish.tree <- read.tree("data/04_fish_tree.tree")
fish.tree  #This tree has 7822 tips and 7821 nodes.

###  Imagine we want to study butterflyfishes. How would we get the butterflyfish clade out of this big tree? It would be cumbersome to drop some 7,000 tips, so it's easiest to use extract.clade. But what node belongs to the butterflyfish family?

### use the function findMRCA in phytools!!

findMRCA(fish.tree, tips=c("Chaetodon_melapterus", "Johnrandallia_nigrirostris"))  

## that is the node you want to use in extract.clade.

tree <- extract.clade(fish.tree, 9164)
tree # this tree has 94 species. 

### Now this is something more easily plotted.
plot(tree)
plot(tree, show.tip.label=FALSE)
axisPhylo() # Shows the scale of the branch lengths.

# So we know that the butterflyfish family has a crown age of about 50 million years.

# Let's say we wanted that crown age to be more precise for downstream analyses. You can get the age using the nodeheights() function.

nodeheight(tree, node=1) ## notice the counting started at 1 instead of #tips+1. 

## Now we want to save this tree for future analyses. Use the "write.tree" function:

write.tree(tree, file="butterflyfishes.tree")


###### Finally, we will do a PGLS. PGLS is a method used to account for phylogenetic relatedness of species when you're interested in correlations between numeric variables. 

## These data for barbets (a king of songbird) include various factors like wing length and song frequency

#link for data: http://www.phytools.org/Cali2017/data/Barbetdata.csv
#link for tree: http://www.phytools.org/Cali2017/data/BarbetTree.nex

data<- read.csv("data/04_Barbetdata.csv")
tree<- read.nexus("data/04_BarbetTree.nexus")

## This tree is in the NEXUS format. The previous trees were in Newick format. You can see the difference if you open the files in a text editor.

## the name.check function ensures that all species in the data are also in the tree and vice-versa.

check <- name.check(tree, data, data.names=data$Species)
check

## There are 9 species in the phylogeny for which we don't have data. Let's just cut these out.

new.tree <- drop.tip(tree, check$tree_not_data)
name.check(new.tree, data, data.names=data$Species) # Now we should be good

## We are going to use the caper package. This works by combining the data and tree into one comparative data object:

comp.data<-comparative.data(new.tree, data, 
                            names.col="Species", vcv.dim=2, 
                            warn.dropped=TRUE)

comp.data ## look at the comp.data object

## We will test if a species' altitude is related to the frequency of its song, while accounting for non-independence due to phylogenetic relatedness.

model1 <- pgls(Lnote~Lnalt, lambda="ML", data=comp.data)
summary(model1)

## Look at the lambda parameter. When lambda=0, this means there is no effect of phylogeny. When it is =1, the evolution of the residual error follows Brownian motion

## You can also plot the liklihood surface of lambda

lm.lk<-pgls.profile(model1, which="lambda")
plot(lm.lk)


