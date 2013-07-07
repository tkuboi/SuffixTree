// SuffixTree construction with Ukkonen's algorithm.

package bioGUI;

import java.io.*;
import java.util.ArrayList;
import java.io.IOException;
import java.lang.Math;

public class SuffixTree {
   private Node bottom, root, prevNode, currNode;
   private ArrayList<String> m_sequences;
   private ArrayList<Node> m_leaves;
   private String m_suffix;
   private Position currEndPos;
   private ActivePoint activepoint;
   private int remainder;
   private boolean fromRoot;
   private int currPos;
   private float totalA;
   private float totalC;
   private float totalG;
   private float totalT;
   private int totalNucleotides;
   private static final String alphabet = "ATCGN$";
   
   private class ActivePoint {
      protected Node node;
      protected char edge;
      protected int length;

      public ActivePoint() {
         this.node = null;
         this.edge = '\0';
         this.length = 0;
      }

      public ActivePoint(Node n, char c, int l) {
         this.node = n;
         this.edge = c;
         this.length = l;      
      }

      public void setNode(Node n) {
         this.node = n;
      }

      public Node getNode() {
         return this.node;
      }

      public void setEdge(char c) {
         this.edge = c;
      }

      public char getEdge() {
         return this.edge;
      }

      public void setLength(int l) {
         this.length = l;
      }

      public int getLength() {
         return this.length;
      }

      public int incrementLength() {
         return this.length++;
      }

      public int decrementLength() {
         return this.length--;
      }
   }

   public SuffixTree() {
      root = null;
      m_sequences = new ArrayList<String>();
      m_leaves = new ArrayList<Node>();
      m_suffix = "";
   }

   public SuffixTree(String str) {
      root = null;
      m_sequences = new ArrayList<String>();
      m_sequences.add(str + "$");
      m_leaves = new ArrayList<Node>();
      m_suffix = "";
   }

   public void addString(String str) {
      m_sequences.add(str + "$");
   }

   public String substr(Position start, Position end) {
      if (start == null)
          System.out.println("start is null!");
      if (end == null)
          System.out.println("end is null!");

      System.out.println(Integer.toString(start.idx) + " : " + Integer.toString(start.position) + " : " + 
                         Integer.toString(end.position) + " size=" + Integer.toString(m_sequences.get(start.idx).length()));
      return m_sequences.get(start.idx).substring(start.position, end.position + 1); // end is inclusive
   }

   public Node buildTree() {
      root = new Node(Node.TYPE_ROOT);
      char c;
      activepoint = new ActivePoint(); 

      int idx = 0;
      for (idx = 0; idx < m_sequences.size(); idx++) {
         int i = 0;
         activepoint.setNode(root);      //set active point to root.
         activepoint.setEdge('\0');
         activepoint.setLength(0);
         m_suffix = "";
         prevNode = root;
         currEndPos = new Position(idx,i); //i will grow to the end of the string.
         remainder = 0;                  //variable to keep track of # of suffix nodes to add.
         fromRoot = true;         //whether or not the active node is originally at root before canonizing.
         currPos = 0;
         while(m_sequences.get(idx).charAt(i) != '$') {
            c = m_sequences.get(idx).charAt(i);
            m_suffix += c;
            System.out.println(i + ": Curr Char is : " + c);
            currEndPos.position = i;
            remainder++;
            update(idx, i, c); //update the tree
            canonize(idx,i); //if the active point is at the end of the edge.
            if (remainder == 0) {
               m_suffix = "";
               activepoint.setNode(root);
               activepoint.setEdge('\0');
               activepoint.setLength(0);
            }
            i++;
         }
         c = m_sequences.get(idx).charAt(i);
         m_suffix += c;
         System.out.println(c);
         currEndPos.position = i;
         System.out.println("----------------------");
         remainder++;
         update(idx, i, c);
         System.out.println("----------------------");
      }
      return root;
   }

   private void update(int idx, int i, char c) {
      currNode = null; //newly created node
      prevNode = root;
      System.out.println("in update");
      boolean endPoint = testAndSplit(idx,i,c);
      while(!endPoint) {
         //System.out.println("prevNode = " + prevNode.getType());
         if (prevNode != null && currNode != null && currNode.getType() == Node.TYPE_INTERNAL) {
            if (prevNode != root && prevNode != currNode) {
               prevNode.setSuffixLink(currNode);
               System.out.println("suffix link created");
            }
            prevNode = currNode;
         }                  
         
         canonize(idx, i);
         endPoint = testAndSplit(idx,i,c);
      }
      if (prevNode != null  && currNode != null && prevNode != root && 
          currNode.getType() == Node.TYPE_INTERNAL && prevNode != currNode) {
         prevNode.setSuffixLink(currNode);
         System.out.println("suffix link created");
      }
   }

   private boolean testAndSplit(int idx, int i, char c) {
      char nextEdge = '\0';
      System.out.println("in testAndSplit");
      if (remainder == 1){
         if (activepoint.getNode().getChildAt(alphabet.indexOf(c)) == null) { 
            //the node does not exist
            currNode = createLeafNode(idx, i, c);
            remainder--;
            activepoint.setEdge(nextEdge); //reset the edge
            return true;
         }
         // if the char is the end of string and the char to be inserted is only $
         if (c == '$') {
            //System.out.println("in 7");
            activepoint.getNode().getChildAt(alphabet.indexOf(c)).addPosition(new Position(idx, i - remainder + 1));
            propagatePositionUp(activepoint.getNode(), new Position(idx, i - remainder + 1));
            remainder--;
            updateActivePoint(c);
            return true;
         }
         else {  //the node does exist
            activepoint.setEdge(c);
            //System.out.println("edge = " + activepoint.getEdge());
            activepoint.incrementLength();
            return true;
         }
      }
      else {  //need to add more than one suffix to the tree
         char e = activepoint.getEdge();
         Node eNode = (e == '\0' ? null : activepoint.getNode().getChildAt(alphabet.indexOf(e)));
         int sIdx = (eNode == null ? idx : eNode.getStart().idx);
         int sPosition = (eNode == null ? i : eNode.getStart().position);

         System.out.println("in else " + activepoint.getEdge());
         System.out.println("edge = " + e + " sPosition = " + sPosition + " length = " + activepoint.getLength() + " remainder = " + remainder);
         // if active_edge is null and active_node does not have the char as its child
         if (e == '\0' && activepoint.getNode().getChildAt(alphabet.indexOf(c)) == null) {
            System.out.println("in 1");
            currNode = createLeafNode(idx, i, c);
            remainder--;
            updateActivePoint(c);
            return false;
         }
         // if active_edge is not null and active_node doe not have the char as its child
         if (e != '\0' && eNode == null) {
            System.out.println("in 2");
            currNode = createLeafNode(idx, i, c);
            remainder--;
            updateActivePoint(c);
            return false;
         }
         // if active_edge is not null and root does not have the char as its child
         if (e != '\0' && root.getChildAt(alphabet.indexOf(c)) == null) {
            System.out.println("in 3");
            currNode = insertInternalNode(idx, i, c);
            remainder--;
            updateActivePoint(c);
            return false;
         }
         // if active_edge is null and root does not have the char as its child
         if (e == '\0' && root.getChildAt(alphabet.indexOf(c)) == null) {
            System.out.println("in 4");
            updateActivePoint(c);
            return false;
         }
         // if the chars do not match
         if (m_sequences.get(sIdx).charAt(sPosition + activepoint.getLength()) != c) {
            System.out.println("in 5");
            currNode = insertInternalNode(idx, i, c);
            remainder--;
            updateActivePoint(c);
            return false;
         }
         // if the char is the end of string and needs to insert char other than $
         if (c == '$' && activepoint.getLength() > 0) {
            System.out.println("in 6");
            currNode = insertInternalNode(idx, i, c);
            remainder--;
            updateActivePoint(c);
            return false;
         }
         // if the char is the end of string and the char to be inserted is only $
         if (c == '$' && activepoint.getLength() == 0) {
            System.out.println("in 7");
            activepoint.getNode().getChildAt(alphabet.indexOf(c)).addPosition(new Position(idx, i - remainder + 1));
            propagatePositionUp(activepoint.getNode(), new Position(idx, i - remainder + 1));
            remainder--;
            updateActivePoint(c);
            return false;
         }
         else {
            System.out.println((sPosition + activepoint.getLength()) + " : " + m_sequences.get(sIdx).charAt(sPosition + activepoint.getLength()) + " sPosition: " + sPosition + " active_length: " + activepoint.getLength());
            
            if (activepoint.getEdge() == '\0')
               activepoint.setEdge(c);
            //increment active_edge
            activepoint.incrementLength();
            return true;
         }
         // endif   
      }
   }

   public Node insertInternalNode(int idx, int i, char c) {
      char e = activepoint.getEdge();
      Node eNode = (e == '\0' ? null : activepoint.getNode().getChildAt(alphabet.indexOf(e)));
      int sIdx = (eNode == null ? idx : eNode.getStart().idx);
      int sPosition = (eNode == null ? i : eNode.getStart().position);
      int ePosition = sPosition + activepoint.getLength() - 1;
      
      //make the internal node explicit
      Node r = new Node(Node.TYPE_INTERNAL,
                new Position(sIdx, sPosition),
                new Position(sIdx, ePosition));
      System.out.println("Active Node is TYPE = " + activepoint.getNode().getType());
      System.out.println("Internal node : edge = " + e + " start = " + sPosition + " end = " + ePosition + " remainder= " + remainder + " length= " + activepoint.getLength());
            
      if (eNode == null)
         System.out.println("eNode is null");
 
      //adjust start index of the original leaf node
      eNode.resetStart(activepoint.getLength());
      System.out.println("eNode begins at : " + eNode.getStart().position);
      //System.out.println("eNode begins at : " + m_sequences.get(sIdx).charAt(eNode.getStart().position));
            
      //add original leaf as a child of the new node
      r.addChildAt(alphabet.indexOf(m_sequences.get(sIdx).charAt(eNode.getStart().position)), eNode);

      if (c == '$' && r.getChildAt(alphabet.indexOf(c)) != null) {
         r.getChildAt(alphabet.indexOf(c)).addPosition(new Position(idx, i - remainder + 1));
      }
      else {
         //create new leaf node
         System.out.println(i + " Create Leaf node : " + c + " position= " + idx + ":" + (i - remainder + 1));
         Node leaf = new Node(Node.TYPE_LEAF, new Position(idx,i), currEndPos, new Position(idx, i - remainder + 1));
         m_leaves.add(leaf);
         //add new leaf as a child of the new node
         r.addChildAt(alphabet.indexOf(c),leaf);
         r.addPosition(new Position(idx, i - remainder + 1)); //for new leaf
         //set parent
         leaf.setParent(r);
      }

      eNode.setParent(r);
      r.setParent(activepoint.getNode());
      //add positions to the new internal node
      addPositions(r, eNode);
      //add new internal node as a child of the parent of the active node
      activepoint.getNode().addChildAt(alphabet.indexOf(activepoint.getEdge()),r);
      propagatePositionUp(activepoint.getNode(), new Position(idx, i - remainder + 1));
      return r;
   }

   public Node createLeafNode(int idx, int i, char c) {
      Node leaf;
      System.out.println(i + " Create Leaf node : " + c + " position= " + idx + ":" + (i - remainder + 1));
      leaf = new Node(Node.TYPE_LEAF,new Position(idx,i),currEndPos,new Position(idx,i - remainder + 1));
      leaf.setParent(activepoint.getNode());
      m_leaves.add(leaf);
      activepoint.getNode().addChildAt(alphabet.indexOf(c),leaf);
      propagatePositionUp(activepoint.getNode(), new Position(idx,i - remainder + 1));
      return activepoint.getNode();
   }

   public void updateActivePoint(char c){
      char nextEdge = '\0';
      //      if the active point is at root
      if (activepoint.getNode() == root) {
         if (m_suffix.length() > 1) {
            m_suffix = m_suffix.substring(1);
            nextEdge = m_suffix.charAt(0);
         }
         else {
            m_suffix = "";
         }
         //System.out.println("0: m_suffix = " + m_suffix);
         //         move the acive_node following rule#1
         activepoint.setEdge(nextEdge);           
         if (m_suffix.length() > 0) {
            activepoint.setLength(m_suffix.length() - 1);
         }
         else {
            activepoint.setLength(0);
         }
         fromRoot = true;
         currPos = 0;
      }
         //      else
      else { // active node is not root
         //         move the acive_node following rule#3
         System.out.println("1: m_suffix = " + m_suffix);
         if (activepoint.getNode().getSuffixLink() != null) {
               System.out.println("have suffix link " + activepoint.getNode().getSuffixLink().getPositionAt(0));
               activepoint.setNode(activepoint.getNode().getSuffixLink());
               fromRoot = false;               
         }
         else {
            activepoint.setNode(root);
            if (fromRoot) {
               if (m_suffix.length() > 1) {
                  m_suffix = m_suffix.substring(1);
               }
               else {
                  m_suffix = "";
               }
            }
            currPos = 0;
            System.out.println("2: m_suffix = " + m_suffix);
         }
         if (!fromRoot && currPos > 0) {
            m_suffix = m_suffix.substring(currPos);
            currPos = 0;
            System.out.println("3: m_suffix = " + m_suffix);
         }
         if (m_suffix.length() > 0) {
            nextEdge = m_suffix.charAt(0);
            activepoint.setLength(m_suffix.length() - 1);
         }
         else {
            activepoint.setLength(0);
         }
         activepoint.setEdge(nextEdge);
         if (activepoint.getNode() == root)
            fromRoot = true;
      }
      System.out.println("next node is : " + activepoint.getNode().getType() + " edge = " + activepoint.getEdge() + " length = " + activepoint.getLength());
   }

   //reset the active edge if activepoint is at the end of the edge
   private void canonize(int idx, int i) {
      Node eNode;
      Node nextNode;
      System.out.println("canonize active_length is : " + activepoint.getLength()); 
      if (activepoint.getEdge() == '\0')
         return;
      eNode = activepoint.getNode().getChildAt(alphabet.indexOf(activepoint.getEdge()));
      int difference = 0;
      //if activepoint at the end of the edge
      while (eNode != null && eNode.getType() == Node.TYPE_INTERNAL && eNode.getLength() <= activepoint.getLength()) {
         //move activepoint to the next node downstream
         nextNode = eNode;
         currPos += eNode.getLength();
         //System.out.println("m_suffix = " + m_suffix);
         System.out.println("eNode length = " + eNode.getLength());
         difference = activepoint.getLength() - eNode.getLength();
         if (difference == 0)
            activepoint.setEdge('\0');
         else
            activepoint.setEdge(m_suffix.charAt(currPos));
         activepoint.setNode(nextNode);
         activepoint.setLength(difference);
         if (activepoint.getEdge() == '\0')
            eNode = null;
         else
            eNode = activepoint.getNode().getChildAt(alphabet.indexOf(activepoint.getEdge()));
         System.out.println("canonize active_node is : " + activepoint.getNode().getType() + " edge = " + activepoint.getEdge() + " length = " + activepoint.getLength());
         System.out.println("m_suffix = " + m_suffix);
      }
   }

   private void addPositions(Node target, Node source) {
      for (Position p : source.getPositions()) {
         target.addPosition(p);
      }
   }

   // add position to all nodes upward in the tree
   public void propagatePositionUp(Node n, Position p) {
      Node currNode = n;
      while(currNode != null && currNode.getType() != Node.TYPE_ROOT) {
         currNode.addPosition(p);
         currNode = currNode.getParent();
      }
   }

   public String getLabel(Node node) {
      int i = 0;
      int idx = node.getPositionAt(i).idx;
      while(i < node.getNumPositions() && idx != node.getStart().idx) {
         i++;
         idx = node.getPositionAt(i).idx;
      }
      return substr(node.getPositionAt(i), node.getEnd());
   }

   public boolean isVisited(ArrayList<Long> list, Node node) {
      long id = node.getId();
      for(long l : list) {
         if(l == id)
            return true;
      }
      return false;
   }

   public static void main(String[] args){
      String str;
      try{
         str = ReadGFF.readFasta("fosmid1.rtf");
         //System.out.println(str);
         //str = "ACTCGGATTATCCAAA";
         //str = "TTTACTTTACGG";
         //str = "GCTCTCGGCTGATGTTTGCTTCTTTTGCTTGTTAAGCTTTGTGCCTACTTTGATTTTCCTGTGCTCCGTGTGTCTCCAGACATTTGTCGTTGGCTCTTAAGGTTTCTGCCTGCTGCCTAGTTTTTCGTTCGTCGTTGATTTCTGTTACCCCTCAAATGCACTTTCAATTGCAGACAGCTTGACGTATATATATGCAGCCGGGTTAGAGCACTATATCTCTCATTTATCTATCTATATACGAACCTACCTATATTTCGTCTGTCGCTGGCTGCTGCCTGCTGGCTGTTGGCTGCTGGCTGGTGCCTCACTTGCACCGAGTCGGAATAGGAATTGGAATTGGAACTGAAATCGGAATCGGAATCCGAGGCACTGGGCTCTGGCTCTGAGCCCCGGCTGCCCACAGGCACTTTGGCTTGTGGGCCGTGGCTCCAACTCCAGATTCGGCTCCAGCTTCCCCAAAGTTTCAGCTGCACCTCGGCTTCTGTCCTGCCCCGGACTTGGACTTGGAGTCGTCGTCAAGTAGCTGCCATGTAAAACGCTTCCTACTTCGGCCTACCAGCGCATTGTAATAATTTGGAAACCTATTGCAGCTCTCGGCCACTACGCAAATTATGAGGTAGCTGGATTTCGTGAATGGATTTGATCATTTGAAATAGTTGAAAACGTCCAATGCCTTATGCAGAAATTGCATCTTCGGTTAACAAAATCTACAATGTCTTCAATTTCGGTTGGCCCATGGTAGTAATATATAGTAGTAATACATAGTTACCTCAATTAAGTTTTCAATTGCCTAACGCGGTGAGCAGAGCATTTCTGCTTCGTTTTCCACATGAACTACTCTAAGTTTTCCTCTCCCATTTCCATTCGAAATGCTTTTGCTGACTGTAATAATACTTTTGTGATGTATTTACTTTAATTAAATTACGTAAAGTGTTACGTGGCCTAAGCACGCACATTTTGAAGACTCCACTGGACGCCCTTCCTCCTCCGAAGTGCCCATATTATTCCAGTTCGCCATACATGTAAGCCACTCCACTTGTGCTTTCATTAGAAAGCCAAAATATTTTTAAGGCTTAAGTCTTCTTGTTGTTTAGCAAGAAGAAAACGAAAACTAATAAATAAATTACACATAAATTCCAAAATTGCCGAAAATTGCCAAGTCTGCCTTAACGGTTAGTCGGCTCAATTGAAAAGCTCTTTGAAAAGCGAGTTTCATTGTCAACGAAATTTATCACCGACTAATGTGAAATTTCACGAATAAGTTTTCAGCGCTGTTGATGGCTTTATTTGCATCGCCAGAAATGATTTAAGCCGCTTCAGCGAAAACTTAATCAAAGACCTAACCCACGAACCTGGACGTTCTCACATAGAAAATGAAAACGAAAGCTGAAAGACAGATGGATCGTTGGATGTGAGCTCTTTTTACACTTGGCTTACACGACGAATGCGGCGTCGGAAAAATTAATATACCAAACAATTACACGGCCCGGCAAAACAATTTTAAGCCCAGTTAATGCCGTCGCACAGTGGCTGCAAGCTACACAAATTCGCTGCGAAAAATTGCCCATAATTGCTAGTGGCCGAGAAATGAAAGCGCTCAAATTTGCATACAAGCAAATGAAGCGCAGCAAAAGGAAAATTGCCACACATACTAAGTGGAAAAGAAATGCAAAACCAATGCGAATGAAAATAATAACAGGGCGCGGTTCGTGAAACTGTCTAATTATTATTTGACATTGGCAATGACATCGAAGCCGAGAGTGCCGCATTTGGATTTTGATTTCAAAAACAGTTGGACTGATGCTGATTAAAAAATATGATGATTATTGGAATAATAATAAACAATGATGTTGGACACAAATGATTACGTCCAGCTTATTGCCTAGAGCATTTCAACAAAGACACTACAATTATTCTAGTGTCTTTGATTTCAAGTTAAAAACATGTATTTTATTTCTAGCTAATTTCAATTCCTACCTAAAAATAGTGAAATATACTTTGATTTGTCATTTCAGCTACTCGTACCACTGTACATTGGCTACAGTCTGCTGAAGCGCACAGTGTTTGAAAGTTAAGACTTTTGCTTTTCTTCCCCGCGACTTTGTTTTTGCTTTTAGAAAGACCAAACAAATGGATAAACAAGAGTTCAAGAGGTTGCCAGCCCCGCCCCTTTGCGCACATGAGGCAAGAAAGAAAAGAAAGTTATTACAAAATTAACAGAGGCTGACAGCTGACAAGGGCAAAGGGATGGGGCAGAAATGGAACTGTACCCAGCTTGAGTGGTGTACGAGAGCCCCGGGCAGAAATTGTGACAAATTTGGGTAAGAAATAGCTGGTAGGAGCGGGACTTACCACCTCCGTTGTTGGCACTGGAGCAGTGCTTGGAGCTGAGGCCCTGACCCGTGTTGGATGCAGAGTTGGGCAGACCGCAGCCACTCATCACGCCGTACTTGATGGGTTTGAAGTGGAAGAAGGAGGGCGAGTAGCCGGATCCGGCGCGGCTCAGGGAGCTGAACTGCTCCTCGCGATCCTCCACGTAGAGGCCCTCCTTGCCCATGGCCGCCAGTTGCCTGCCCTTGGGCATGAACATCACCAGGAAGACAGTGGCGGAAGTGGCTACCAGGCCGAAGGCCACGCAAGCGTCCTTGTGCCGCTCGGCCACGGCCAGTCCGCATAGCATCCATCCCAGCCAGATGGGGATGGCTCCTCCGATGGCCAGCCCGATGTAAGTGGCCTCACGGTAGTTGTCCCGGATGCCGCGCGACTTAATGGCCAGCACGGCGATGAAGACGATCAGAAAGATAATGTAGATCAGCGAGAAGAGCAGCTCCGAGAACTGGGTCTTGCAGAGCGGTATTAGCACGGTGCTGACCGCCGCTATCCGGGTGTAGATCTCCGGAGTGCCGTCCAGCGTGGTGTAGGACGTGGGGTAAAACAGTGCCGAGTAGTTGGTCTGCGAGGCCACGGTGGTGCTAAGAAAGCCACTGCCCATCACAGGCACACTGGTGGTGTACACCTCCGGTGGTTGAGTGAGCAGCCACTGACCTCCGATCGCGACCTGGATGAGCAGCGCAAAGAGGAGCAGCAGTCCCTGGTATGGAGCCGGCAGATACACCCCTCCATTCAGACTAATCAGAAATACGCACTTGACCAGCAGGGCGGCAAACACCAGGGCGTAGGCCACGCCCACTCCGAAGCGGATGGCCCCGCAGCTGATCAACGAGGGCTGAGCGGTGATGATGGCCCCTAAACTAGCGCAGGCGAAAAGGCCGAGGAGCAGCATCTGGCCCAGAAAAAGGTGCCTTCGCGACGGCGACGTCCGCCAGGCCTTGAAGAGCACAAAGATCTCAAAGGCTGCCATCATAAGCATGGTGAGGGTGGCCAGGACTAGCACCGGCACCACCCACGGCTCCTTGCGTAGGCCTGCGAAGGGCAAGTGGGTGGACACTATGCTGGGGGGCGTGCTGTTGGAGCTGCTCCCGGGAGTCACCCGGCCAGTGGTGCCGGCTGCTCGCGCCTTCGCCCCCGAACTCCGACGAGGCGTGGCTGGTGTGGTGGTGGACGTGGTGCTGCTAGAACTGCTGGCCGACAGTGTGGGTACCACAAAGAATTCCGTCGACATTTCGAGACGGGCGGCACTGGAGACATTGCCACCCGGAAAATTGAGATCCAATTGCTTGGCAGTTCGTTGGGACTGAGGAGCTGCCGGCGTAGTGGTACTAGTGGTGGTACTTGTGGTACTCGTAGTAGTAGTAGTGGTGCTGGTGGTATTACGGCTACTCGGCTTGAGCACCTGGTAGGCGTTACGATCGGTCTCCTGGCTAGTCACCAGGCGCACCCGCCGCCAAATGTTGTTTCCGTTAACGCTGAGGCTGCGCGACGCCGGCGTGGTGGCGATCAGCTGCAGGATCAGCAGCAGGATCAGCAGCAGGTGGCGTCTACCGGCCATTGAAAGTGCTGACCTTCGATGCCCTGGATCTCTGATCAAGTAATCAAACGCGGCACATCTCTCGCGAACGCACACTAACCACGCACACGCGATCGCGATTCTTCAGGCTCGTGAATTCTTCAGTCTGAGCAGGCCCTCGATCTTCAGTTCTGCGGTCTTATTCTTCTCGATTCTCGCGCTAAAACGTCCGCTTGCTCCGATCGCCCAACTGCAACTGACTGGCGGCTCTGCAAGTTGGGCTTTGGACCGGCTTAGTCGTTTGAGACCAGTTTGACGGCAGCCCAGTCGAACTAACACTTGTTTGAGTTCGATGTGCGTTTGCTAGCCCGGCTAACACTTTCAGGAACGCAATATTATAGCGTTGAGGCGTGGGAGGCAGTTGGGCTGTAAGGGGTTTATATGATTTTTGGCAACTTAATATTCCCCGAGCATCTAAATACTTATTAAATGTAGTTACTCGATCTTATGGGTTCTAGGGCTGAAGCATAGGTTCAACAACTACTGTTCGCAATCGCTCGACTACTCGGCGATATCTATGGATTCGATCTCTAATATCTATATCTATATATTATATCTATAATTTGATATTTGTCTTTAATTCCCGTACTATCCGTTGGTGTAAAAACTTTGGGAAGAGAAACCCATATAATCCAGGAACTCATTTACGTATTTATGAAGTTAGAAATAGTTAACAAATAAAATACATGTTCCATTTGATAATAATTCCAATTTGATTTCAGATTATGGGATACCATCGGTCTACCACCGCCTTTTCCCTGCGCAGTTGGTTGCAGCCATCCACAAGCAGACGCGATCCACAAACTTAGCCCAGGATCCACACTCAAACCGATTGCGGGATGGAGAACACTTTCATTCGTTCTGCATATGCGGAATTTCAGTTTATGGGCTCGCTGCTCGACTTGGCTGGTAGGAAAGTGCTACGGGCCATTATCTATTAGCCAACCGCTTTTGCGGTCTCCTCGCGGAGCCGCTCTCCTCGAAAACGACCCAACAAAAGGCCAAATAGCGAAATCCGCGCACAAAAGGCGATTGTAAATCAAAAGAGCGCTAATGATATATCGATATACAAAGTGTTGAGGCGGGAAAGGCAACGCCAATGGCGTTAGAAACGGCAATAGCAGTAGAAACTCGCTCTGATTCCTTAGGTATTATTGTGTTCGACGCGCGCTTTCAGCTAACCACAGTTAGAATTCAGATCACCTGAATGAGTGGGCAAAAAGGCGTTAGAAAAGGTACATATTGCATAAAATCGGCGTTCAAGTTGGTCTAAGCCGTAATCACACATGGTCTTTAGGGACTCCCCACAGCTACTGACCTAAATAAATATCTCACAGCTATTTCCTATATCTTCGAGCTACGTATATAACAATAGATAATAATGATGTTATAGGAATCTCGAAGGTATGTAATCTTTGAACTTATGTAAGATGTACTTATATGGCAATCCGATTTAACTAGCTCTGATGCCAAATTGTCCTGTGAGCGGTAAAGGCCAACTCTCTGAAGAATCTTGAGGTACAACCACTGGTCGAGTACGAGTCTGCGTTTTGGCCCAAGCGGACCGTATTGGCCACTAGACATTTATCTCCATATTTGGCGCGACCCGCAATACGCAGCCATTAGCAGCAACAAATGTTATGTAGATCGTGTCGACGGCGACTGGGCTCCACCGCCACCCCCACCCCCTTAGGACCGGACTAACTGGCCACCGCCCGACTTGCCAACTTACGGTATCGGTAGTGACGATCTGCGAGGGTCTGGAGGGTCGGAGTCCGGCATAGCGTTCTGTATTGACATGGAGGCGGGTAAATAACCTTTCGGCATGGCCCGAATACAGAGCAGTCAAATTTGTGACTGGCCATAACCGACCATAACTCAACTCGGCTTCTGCTGCGATCATAATCTAGAGCGATACGGGGGTATATACTATATGCGTGCCAGTCGCTCTTCCCAGTGCCTGTTGGCTAAGTAGATCTTATATTAAATAATATTGTATTACGAAGCGAGCAGGGATGATGAGGCAGCGTAATTTACGATCTTTTCAATTATTACTTGTCGCCGCTGAATAAAATTCCCTTCACGGTGTGCTTAATCCTGGCATTAGCTTATGCCGCTGTAGGCCGATTTCCATGATCCGTTTAATTAAAGCTGTCACTTGCCTGTTAGCCTTGCCAAACGGCCAGCTTCCCTTGGCAACAGTCCCAGTAACAGTGCTTTATGTAATCCCGCGGTCTTAGAGCCATCTTAGCTAATTAAAAGGCATTCTAATATTACGTTTCGCCTTTCTCCCGACCGCCAGATCACGTTTCACGACATACCCGAGATCCGAGACACATTGCGTTAATAACCCCTGACTGAGCACAAGATCAAGGCAGCGCTGTGGCTCGTGGTCAGGCCAAAAGTCAATAGGCCCACACCTTGGCCAAAAAATTGGTCCTGCTCTCGAAATAACCAAGTTTCATGTTGACTCATACCTGGCAGCTAGTACTATTGATTCCTACTGACTCGTATATTGGAAAAATTGAAACAATACCAAAACAAGAAAGAAAGCTAGCTTTGGCTAGCCAACGCTATACATGTATATAAAAACACTTTATGATCAAAAACAACGAAACTACTTTTTCTTTTTATTATTTTCCATTAATTTTCCGAATGTTCCTATGGCAGCTTTATGATATAATTGTCCGATTTTGGGAAAGTCGCATTTAAATAACTTTTAAGTTGGGAGACTATCTTGCGCATCAACGGACGAGCGGCTGTTAATACTGATCAAGAATATAAATAAGCAAAAGTCGCTTGCACTGCACTGCAAGCTGACTGAAATATCAATATGAGGAGCATAGAAAGGAACTTAAGTAAAATTCATGTATATTCAGAATCCCTGGATTCTGAGCATAAATGATCGTTATCGTTTTCCGTAGTGCTCTCATTACACGACTTGTGTATCAACTTTCTTTATGCCTAAAAATTATTCAAATTACTCGATGAGCCTAGAATATTGGCACGACTTCGCAGGGAGAAGTAACCAAACTGGCCGTGACCTCGCCGCCTAGATGGGCCAATACTACTTTGGTCGAATAGACCGCCGGCGCCAAAGATTGCGAATCCGTTCCAAAATTGATTGATTCCGGGTGATTTGAATGTCTGGCCCGAGACCTGGACAACCTCCGGGACCGGGACCTGGAACCCCGAGAACTTCTAGATTCTAGAGGATGAGGAAGTGCCGCCAGTTGGCTTGGCCAACGAAATGGGTCCGATAATCGGGCCGGACCTCTGCTTGCGGCTCTTTTTAGGCGCCATAATTTATATGCCAGAAGACAAGACGGTGCCGTCTGGGTCAAGTGTCAATCGGCTTACGATTTGCAAATGTTTCCGCTCATAATTCTATGGTGCAATTAATTTTTGGACACATCGATCGTTGAATTTTTGAGCAAGCAAAGACGCCGCCCGTAAATCTTCGGGGATGCATATTTAATTTAAATTGAAACCAGTTCTGCCAGTGCATTGAACCTGCCACAATCGGAATGAGCGGGACTCGGATGCCAAATATTCAGAGACAAAAGCCTTTTGGCGATCTTATTTATTTATCAATATATAGTGCGGCTCGCAACCCTTTTCCCCATTTCGTTCGGCAACATTTTCCCCGGGTTCATTGTCGTCCCTAAAAGAGCGAATTCTGATTTATAAGCGGGTCTTGGCCAAATGGCCGTCTTCTCCCTCATCCCCCCACCAGAACCCCTTCTTCGCGGCGGGGCTGCGCTTTAATCACGTGCGTGCAATGCGCATATGATGATGCGCCCCAATCCGAATCCCCCTGCTTGAGCCTTCAGTATGGCTCGATCCCTGGTCCCCGATTCCATCCCATCCCATCCCATACCATCCCAACCCAACCCAGAAGCGATCCCGACCCCTCCGATGCTCTTGACGATGCTGGGATGCCGCCCTGGACTTGCTTTTGTTACTGGTGTTGCAAGTGCTGGCAGGCCAAACACCCAGTGAACATTGCTCAACGGAGGAGTTCGCTTTCTTGCAACTCACAGTGGTGGGCAGGTGCTTAGGTACTCGTTTTCATACGCAAGCAAAATAATAAAGATGGTGAATGTGGAATATAGAAATGATAAATGCATTTCTGAGCCACATGTTACATATCTTATTATTTTGATCAGGATGGCCGCCGACATTTCTTATTAAATTATTATTCCCATCCCTGACCAATGAAGGCGAAAACGTAAAGTTCTGTTTTTCCATTATACTCTTGTATATCTACGCAGCATTCCATTGTTTGAGAGAATCATTTTTCGTGACACTTCTCTTCTGGGGCCTAGTTCGCTGACCTTGCATATTATATATCTCCTATTATTTTGATCAGGATGGTCTCCGACATTTCTTATTAAATTATTAAATAATTATTCCCGACCCTGACCAATGAAGGCGAAAACGTAAAGTTGTTTTTCCATTATACTCAATACTTGTATATCTGCGCAGCATTCCATTGTTTGAGAGAATCATTTTTCGTGACATTTCTCTTCTGGGGCCTAGTTCGCTGACCTCCGCTGTGCTTCCCACTTTTGGCTTATTCCATTATGCAGAATGCATATAAATGCAAGTGCAAGTATCGAAGTTAGTACGGGAGTGTTGGCCCGACTTGGCTGCCGTTGTGTGCGAGCACGCATGCGTGAAGCGTCCGTTTTAATCGTTGCCTGGCCATCGAATATTTCTTTGGCCACGACGCTCCGACGGCGATGACGACGTGGAGAGTTTGCTCGTAAATCCGGCGTATAAAAGGCACATAGCTCCCAATGGTGCCCCAATGGTGCTGGCCGGGAGGGGGAGCAGGGCAGGGGGTGTCGATTGTGCTGTGGGGAGTGGTGGGGTCAGCCCCGATGGCGAACTTGGTCGAGCAAGGCATAAGGCCTTCAGATGGCGGCTAACTCGGCTCTAAGGTGCCCGAGCGCAGTGGTCCGCCTGACCAAAGGCTGATTTACTATGATCCGCTTGAGTTTTATGTACAGAGAGATAAAATTTTTAAGATATTGACACCTTGTTAGGCGTTTGAGAGACAGCCGGAGCCGGGATTGCCGAGGGACTGCCATAAGGAAAAGAACTCTTCTTGGCCATAATCTCCACCGGGAACGCCCGCCGATCGCACCTTCCCTGGATTCCCATTCCCATTCGCATTCGCAGTTGCCCTTTTGCGTGGCTAATTGTGGTTAAAATTAAGTGCATTTCTCGCCGTCAAATGCATAAACTTGGGCCAGATAATCAGGACGACAGAGGTTCTCTTGGCCTGTGCGTGGGTGTCAAAGGCGGGCTTTGTTAAATTTATGGAATAATACATATCTTTACGGCTTCATTAGCGAAAACCCGTCCCCGTCCCCGTCTCCACAAGAAAAGCCATTCGACACGCTGTCTGCATTAGCGCCAAAAGGATCGGCATCGACATAGTCATCAAGCGTCATCCTCTTCGTCTTCGTCTCCAAGACTGCGTTGTCTATCTATCTGTCACTACGGAGAGTCGGTCCTCATTCTCGAGGCCCTGAGCCCGCGATTCTTTGTATGCACACTTCAAAGGGGTTCAGCTCGGGTCGCTTCGGTTGAGGAGAATGCTTCAGAGTAAGGAGATTCTGATTCCCAGTTCACGAAATCTGCATGTTTAATATATCTGCTGCTATTACTTGTCTAATAATAATAAGGTTTTAACACACCCAACTAAAATACAATTATGATTTCATTCATAGTATCTGGTATGGTCAAAAAATGCCATTTTAGTAGAACACATTTCAAAAACATTCTAACTTGATTAGATTGCAAAGTAGGTCATTCTTCTATTCAAATCAATTCATTTTGCGTTGAGCACTTTTATAAGTGCACCCCCAATTCGGTCACTTTTTCGGCGGAGAGCAGTCGGGTTGGTTCTTAGCTTCCCAGATTGATTGAGCCTTTTATTAAAGCGCACGTACAGCGTATTAAATAAATTGTTAGTTACGTTCCGGAACGGAACGGAATCGAAGGGAACGAGCCCGTCCCAGTCTCCAGACTCCAAGCCCAAACCCCAACCAAAGCCCTTCTCCTACGTCAGATTTTGTGGCAAGAGTGCTGTGCTGCGGCTCGGTATGGTGATGGTGATGGTGTCATATGGCACAACGCCCCCGAAGGTGTTGTGCGGGAAAATCTGAGGAACATGGCGGGAATGGACACAACGAAAGGAATTCCTTTTCCTACAACGTCCTACAACTTTGTGCGCGGCGAGCCGCGTCTGTGGGCCCGGGTCTGTCCAGATAGATGTCTATGCATTTTGCTTTATGGTCCAGCCACGTATCAATTATGTGCCAAAGTGCTAATGGCTGTTGACGGCTTCCAGCGTCCATTGTGGGTGTTGCCCCAAAGGAACATTCCTCGGACGACGGACGACGGACGACGTACGTATGCCGCAAAGGGTTTCCATCTGCCGCGGCTAAATTAGTTTTCCACTCCTCTGCACTGCGATGCAAGCGATTTCCTTCGGGCAAAAAGTGAAAACCAGTTTTAATTACCTATCCGCTTCGGCGGCGAATATCGCGTGTCGATCTGCAAGCGGAGCAATTCAATTATGACACAAATCCATCGCATGCCAAGCCAAATCCAATTCGGCATTAGATATCCAAAAACCTGCAAAGCCCAGTTAACGCAGTCGTGGCCAACTTTTCAAAGGTGTGTACTTTCTAATAGGTAAATGCTTACACTGGATAACAATAAACCTCTTAATATACTTAACTATTTTTACAAAACCAAAGTATTACCCTCTGTAACCTACTACCAAATGGAATTTAGCAGGCGTTTTGTGGCTTGTCAGTTTTTTAATGGGACATTAATTGAAGACCATAAATATCTAATTGCCCAATTGCTTAGTCTTCGTTGGCTTATATTAAATATGTAAAATGATCCCCGTTCAGGCTACCACCACTCATAGAGTTTCCATTGAACATAGCATAGATACATAGATTCTGGCCAGTATCGAGACAAACATTT";
         str = "GCTCTCGGCTGATGTTTGCTTCTTTTGCTTGTTAA";
         SuffixTree sf = new SuffixTree(str);
         //sf.addString("TTTGGATAATCCGAGT");
         //sf.addString("AAAAAAAAAA");
         //sf.addString("ACGG");
         Node r = sf.buildTree();
         System.out.println(r.getNumChildren());
         for (int i =0; i < 5; i++) {
            if (r.getChildAt(i) != null) {
               System.out.println(i + " : " + sf.substr(r.getChildAt(i).getStart(),r.getChildAt(i).getEnd()));
               System.out.println("length : " + r.getChildAt(i).getLength());
               for(Position p : r.getChildAt(i).getPositions())
                  System.out.println("position : " + p.toString());
               for (Node n : r.getChildAt(i).getChildren()) {
                  if (n != null) {
                     System.out.println(sf.substr(n.getStart(),n.getEnd()));
                     System.out.println("length : " + (r.getChildAt(i).getLength() + n.getLength()));
                     for(Position p : n.getPositions())
                        System.out.println("       position : " + p.toString());
                  }
               }
            }
         }
         System.out.println("number of leaves = " + Integer.toString(sf.m_leaves.size()));
      }
      catch (Exception e){
         System.out.println(e);
      }
   }
}

