package com.example.fluiddynamicsjavafx;

import java.util.List;

public class NeighbourThread extends Thread{
    List<KdTree.Node> particles;
    public NeighbourThread(List<KdTree.Node> particles){
        this.particles = particles;
    }
    public void run(){
        for(KdTree.Node particle : particles){
            Fluid2D.neighbors[particle.id] = (Fluid2D.tree.rangeSearch(particle,Fluid2D.smoothingLength));
        }
    }
}
