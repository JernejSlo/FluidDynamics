package com.example.fluiddynamicsjavafx;

import javafx.scene.paint.Color;

import java.util.*;

public class KdTree {
    private int dimensions_;
    private double radius;
    private Node root_ = null;
    private Node best_ = null;
    private ArrayList<Node> collisions = null;
    private double bestDistance_ = radius*radius;
    private int visited_ = 0;
    ArrayList<Node> neighbours;

    public KdTree(int dimensions, List<Node> nodes) {
        dimensions_ = dimensions;
        root_ = makeTree(nodes, 0, nodes.size(), 0);
    }


    public ArrayList<Node> rangeSearch(Node base, double radius) {
        // Create a list to store the found nodes
        ArrayList<Node> range = new ArrayList<>();

        // Call the recursive rangeSearch helper method
        rangeSearch(root_, base, radius, range);

        return range;
    }

    private void rangeSearch(Node root, Node base, double radius, ArrayList<Node> range) {
        // Base case: if the current node is null, return
        if (root == null) {
            return;
        }

        // Calculate the distance between the current node and the search point
        double d = getDistance(base.coords_[0],base.coords_[1],root.coords_[0],root.coords_[1]);
        if (d <= radius && base.id != root.id) {
            range.add(root);
        }// Recursively search the left and right subtrees, if they exist
        rangeSearch(root.left_, base, radius, range);
        rangeSearch(root.right_, base, radius, range);
    }


        public ArrayList<Node> findCollisions2(Node target, double radius){
        double start = System.currentTimeMillis();
        neighbours = new ArrayList<Node>();
        if (root_ == null)
            throw new IllegalStateException("Tree is empty!");
        Node currPosition = root_;

        int index = 0;
        double dx = currPosition.distance(target);
        while (dx > radius){
            System.out.println("dx " +dx);
            dx = currPosition.distance(target);
            currPosition = currPosition.left_;
            if (currPosition == null){
                return neighbours;
            }
            index = (index + 1) % dimensions_;
        }
        System.out.println("distance "+currPosition.distance(target));
        recAddAllToNeighbours(currPosition, target.id);
        //System.out.println("Found where they collide in " + (System.currentTimeMillis()-start) + "miliseconds");
        return neighbours;
    }

    private void recAddAllToNeighbours(Node curr,int id){
        if (curr.id != id){
        neighbours.add(curr);
        }
        if (curr.left_ != null){
            recAddAllToNeighbours(curr.left_,id);
        }
        if (curr.right_ != null){
            recAddAllToNeighbours(curr.right_,id);
        }
    }

    public void displayNode(Node root){
        if (root == null){
            root = root_;
        }
        System.out.println(root);
        System.out.println("<-  ->");
        System.out.println(root.left_+" | "+root.right_);
    }


    public ArrayList<Node> findCollisions(Node target,double smoothingLength) {
        if (root_ == null)
            throw new IllegalStateException("Tree is empty!");
        collisions = new ArrayList<>();
        visited_ = 0;
        bestDistance_ = radius*radius;
        nearestInRadius(root_, target, 0);

        System.out.println("number of collisions " + collisions.size());

        for (int i = 0; i < collisions.size(); i++) {
            if (target.id == collisions.get(i).id){
                collisions.remove(i);
            }
            else {
                System.out.println("collided with " + collisions.get(i).id);
            }
        }

        return collisions;
    }

    private void nearestInRadius(Node root, Node target, int index) {
        if (root == null)
            return;
        ++visited_;
        double d = root.distance(target);
        System.out.println(d);
        if (d <= bestDistance_) {
            System.out.println("its here");
            collisions.add(root);
            best_ = root;
        }
        nearest(root.right_, target, index);
        nearest(root.left_, target, index);
    }

    public Node findNearest(Node target) {
        if (root_ == null)
            throw new IllegalStateException("Tree is empty!");
        best_ = null;
        visited_ = 0;
        bestDistance_ = 0;
        nearest(root_, target, 0);
        System.out.println(best_);
        return best_;
    }

    private void nearest(Node root, Node target, int index) {
        if (root == null)
            return;
        ++visited_;
        double d = root.distance(target);
        if (best_ == null || d < bestDistance_) {
            bestDistance_ = d;
            best_ = root;
        }
        if (bestDistance_ == 0)
            return;
        double dx = root.get(index) - target.get(index);
        index = (index + 1) % dimensions_;
        nearest(dx > 0 ? root.left_ : root.right_, target, index);
        if (dx * dx >= bestDistance_)
            return;
        nearest(dx > 0 ? root.right_ : root.left_, target, index);
    }

    public int visited() {
        return visited_;
    }
    public double getDistance(double particleX,double particleY,double neighbourX,double neighbourY){
        double dx = particleX - neighbourX;
        double dy = particleY - neighbourY;
        double distance = Math.sqrt(dx * dx + dy * dy);


        return distance;
    }
    public double distance(double[] point) {
        double dist = Math.sqrt(Math.pow(point[0],2 )+Math.pow(point[1],2 ));
        return dist;
    }



    private Node makeTree(List<Node> nodes, int begin, int end, int index) {
        if (end <= begin)
            return null;
        int n = begin + (end - begin)/2;
        Node node = QuickSelect.select(nodes, begin, end - 1, n, new NodeComparator(index));
        index = (index + 1) % dimensions_;
        node.id = n;
        node.left_ = makeTree(nodes, begin, n, index);
        node.right_ = makeTree(nodes, n + 1, end, index);
        return node;
    }

    public void setRadius(double radius) {
        this.radius = radius;
    }

    private static class NodeComparator implements Comparator<Node> {
        private int index_;

        private NodeComparator(int index) {
            index_ = index;
        }
        public int compare(Node n1, Node n2) {
            return Double.compare(n1.get(index_), n2.get(index_));
        }
    }

    public static class Node {
        private double[] coords_;

        public Color getColor() {
            return color;
        }

        public void setColor(Color color) {
            this.color = color;
        }

        private Color color = Color.LIGHTSKYBLUE;
        public int id;
        private int radius = Fluid2D.PARTICLE_RADIUS;
        private double [] velocity;
        private double density;
        private double pressure;
        private Node left_ = null;
        private Node right_ = null;

        public Node(double[] coords) {
            coords_ = coords;
        }
        public Node(int id) {
            this.id = id;
        }
        public Node(double x, double y) {
            this(new double[]{x, y});
        }
        public Node(double x, double y, double z) {
            this(new double[]{x, y, z});
        }
        double get(int index) {
            return coords_[index];
        }
        double distance(Node node) {
            double dist = 0;
            for (int i = 0; i < coords_.length; ++i) {
                double d = coords_[i] - node.coords_[i];
                dist += d * d;
            }
            return dist;
        }
        public String toString() {
            StringBuilder s = new StringBuilder("(");
            for (int i = 0; i < coords_.length; ++i) {
                if (i > 0)
                    s.append(", ");
                s.append(coords_[i]);
            }
            s.append(')');
            return s.toString() + " with id: " + id;
        }
        public void show(){
            System.out.println("For particle " + this.id +" the values are:");
            System.out.println("coordinates: ("+this.coords_[0]+", "+this.coords_[1]+")");
            System.out.println("velocity: ("+this.velocity[0]+", "+this.velocity[1]+")");
            System.out.println("density: "+this.density);
            System.out.println("pressure: "+this.pressure);
        }

        public double[] getVelocity() {
            return velocity;
        }

        public double getDensity() {
            return density;
        }

        public void setDensity(double density) {
            this.density = density;
        }

        public double [] getCoords_() {
            return coords_;
        }

        public void setCoords_(double [] coords_) {
            this.coords_ = coords_;
        }

        public double getPressure() {
            return pressure;
        }

        public void setPressure(double pressure) {
            this.pressure = pressure;
        }

        public void setVelocity(double[] velocity) {
            this.velocity = velocity;
        }

        public int getRadius() {
            return radius;
        }

        public void setRadius(int radius) {
            this.radius = radius;
        }
    }
}
