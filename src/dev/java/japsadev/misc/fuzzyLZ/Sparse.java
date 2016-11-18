/*
 * Copyright (c) David Powell <david@drp.id.au>
 * 
 * 
 * This file is part of FuzzyLZ
 * 
 * FuzzyLZ is a program orginally intended for the compression of DNA sequeces.
 * It can be viewed as a compression model like Lempel-Ziv 77, but instead of
 * exact matches, allowing matches that contain inserts/deletes/mismatches.
 *  
 */

/*
 * This class implements a sparse array to keep track of the active cells in the
 * DPA
 */

package japsadev.misc.fuzzyLZ;

import java.io.*;

import java.util.*;

import japsadev.misc.common.*;
import japsadev.misc.dnaPlatform.gui.MainFrame;

class Sparse implements Serializable {
	static final long serialVersionUID = MainFrame.serialVersionUID;

	class Link implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		int start, end, age;

		@SuppressWarnings("rawtypes")
		ArrayList o;

		Link next;

		Link prev;

		Link(int start, int end) {
			this.start = start;
			this.end = end;
			age = 0;
			next = null;
			prev = null;
			o = new ArrayList<Has_Value>(end - start + 10);
		}

		@SuppressWarnings("unchecked")
		void realloc() {
			o.ensureCapacity(end - start);

			while (o.size() > end - start)
				delObj((Has_Value) o.remove(o.size() - 1));

			while (o.size() < end - start)
				o.add(newObj());
		}

		public String toString() {
			return "l(" + start + "->" + end + " " + age + ") ";
		}
	}

	/**
	 * Class Iterate is used by a number of functions to store state information
	 * about the current position
	 */
	static class Iterate {
		Link curr;

		int i;

		Has_Value o;

		public String toString() {
			return "Link:" + curr + " i=" + i + " o=" + o;
		}
	}

	final Has_Value type;

	Link head; // Keep these sorted by link.start

	Link tail;

	int numLinks;

	int numObj;

	long links_created = 0, newObj = 0;

	LinkedList<Has_Value> unusedObj; // Storage of FSM not being used. So we
										// don't have to

	// keep destroying/creating new ones.

	Sparse(Has_Value type) {
		this.type = type;
		unusedObj = new LinkedList<Has_Value>();
	}

	/*
	 * void clear() { Link l = head; while (l!=null) { l.clear(); Link t =
	 * l.next; l.next = null; l.prev = null; l = t; } }
	 */

	private Link insert(Link ins, Link prev) {
		links_created++;
		if (prev == null) {
			if (head != null)
				head.prev = ins;
			ins.next = head;
			head = ins;
		} else {
			// Insert ins into double-linked list.
			if (prev.next != null)
				prev.next.prev = ins;
			ins.next = prev.next;
			ins.prev = prev;
			prev.next = ins;
		}
		if (ins.next == null)
			tail = ins;

		return ins;
	}

	private Link copy(Link l, Link prev, int start, int end) {
		if (l == null)
			l = insert(new Link(start, end), prev);
		else {
			l.start = start;
			l.end = end;
		}
		return l;
	}

	void copy_cull_cut(Sparse s, int dir, int min_age, double cutOff, int maxGap) {
		Link l = head, prev = null;

		sanity();
		for (Link ol = s.head; ol != null; ol = ol.next) {
			if (ol.age < min_age) {
				// This link is young, copy it completely
				int start = ol.start;
				int end = ol.end;
				if (dir == 1)
					end++;
				else
					start = (start > 0 ? start - 1 : 0);
				l = copy(l, prev, start, end);
				l.age = ol.age;
				prev = l;
				l = l.next;
				continue;
			}

			int start = -1, end = -1;
			for (int i = ol.start; i < ol.end; i++) {
				Has_Value o = (Has_Value) ol.o.get(i - ol.start);
				if (o.get_val() <= cutOff) {
					if (start < 0) {
						// Start of a new significant part
						start = MyMath.max2(i - maxGap, ol.start);
						i = MyMath.min2(i + maxGap, ol.end);
						end = i;
					} else {
						end = i + 1;
					}
					continue;
				}

				// This Obj has a high msgLen

				if (start >= 0 && i - end > maxGap) {
					// We have a significant region with a big gap after it.
					if (dir == 1)
						end++;
					else
						start = (start > 0 ? start - 1 : 0);
					l = copy(l, prev, start, end);
					l.age = ol.age;
					prev = l;
					l = l.next;

					start = end = -1;
				}
			}

			// Any remaining bit at the end?
			if (start >= 0) {
				if (dir == 1)
					end++;
				else
					start = (start > 0 ? start - 1 : 0);
				l = copy(l, prev, start, end);
				l.age = ol.age;
				prev = l;
				l = l.next;
			}
		}

		// Must remove anything left at the end of l (this sparse list).
		while (l != null) {
			Link n = l.next;
			remove(l);
			l = n;
		}
	}

	// This copies the format of sparse array 's' into this sparse array.
	// It is assumed no object data is in the array. And checkAlloc() must be
	// called before the array is used.
	// 'dir' can be -1,0,1 to specify which direction to expand the sparse array
	void duplicateFrom(Sparse s, int dir) {
		Link l;
		Link l2 = head, prev = null;
		for (l = s.head; l != null; l = l.next) {
			int ns = l.start - (dir == -1 && l.start > 0 ? 1 : 0);
			int ne = l.end + (dir == 1 ? 1 : 0);
			// System.err.println("Duplicate link ("+ns+" -> +"+ne+") l2="+l2);
			if (l2 == null) {
				l2 = insert(new Link(ns, ne), prev);
			} else {
				l2.start = ns;
				l2.end = ne;
			}
			l2.age = l.age;

			prev = l2;
			l2 = l2.next;
		}

		// Must remove anything left in l2.
		while (l2 != null) {
			Link n = l2.next;
			remove(l2);
			l2 = n;
		}

		sanity();
	}

	// Join any touching or overlapping sparse elements. Assumes no useful data
	// in the array
	void join() {
		Link l = head, p = null;

		while (l != null) {
			if (p == null) {
				p = l;
				l = l.next;
				continue;
			}

			// Join any overlapping regions
			if (p.end >= l.start) {
				p.start = MyMath.min2(p.start, l.start);
				p.end = MyMath.max2(p.end, l.end);

				remove(l);
				l = p.next;
			} else {
				p = l;
				l = l.next;
			}
		}
		sanity();
		sanity2();
	}

	// Remove a link from the sparse list
	void remove(Link l) {
		if (l.prev == null) {
			Misc.my_assert(head == l, "head!=l in Sparse.remove()");
			head = l.next;
		} else
			l.prev.next = l.next;

		if (l.next == null) {
			Misc.my_assert(tail == l, "tail!=l in Sparse.remove()");
			tail = l.prev;
		} else
			l.next.prev = l.prev;

		// Now set its size to 0, and realloc() so all Object are reclaimed in
		// our temp list
		l.end = l.start;
		l.realloc();
	}

	// Add a new one from 'n1 to n2'
	// NOTE: This function assumes that no useful object data is currently in
	// the array.
	// This function does _not_ allocate room in links.
	// After using this function, checkAlloc() must be called before using this
	// sparse array
	void add(int n1, int n2) {
		if (n1 < 0)
			n1 = 0;
		if (n2 <= n1)
			return;

		Link l = head, p = null;
		while (l != null) {
			if (n2 < l.start) {
				// Put new one before l
				Link nl = insert(new Link(n1, n2), p);
				nl.age = 0;
				return;
			}

			if (n1 < l.end) {
				// Join new one into l
				l.start = MyMath.min2(n1, l.start);
				l.end = MyMath.max2(n2, l.end);
				// l.age=0;
				return;
			}

			p = l;
			l = l.next;
		}

		// Put the new one at the end.
		Link nl = insert(new Link(n1, n2), p);
		nl.age = 0;
	}

	void checkAlloc() {
		sanity2();
		for (Link l = head; l != null; l = l.next) {
			l.realloc();
		}
	}

	// Remove elements from the start and end of each sparse list that are >
	// cutoff
	// Skip any sparse lists that have an age < minage
	// Not used by current implementation.
	void cull(int min_age, double cutoff) {
		Link l = head;
		while (l != null) {
			if (l.age < min_age) {
				l = l.next;
				continue;
			}

			int len = l.end - l.start;
			int s, e;
			for (s = 0; s < len; s++) { // Find first element < cutoff
				Has_Value o = (Has_Value) l.o.get(s);
				if (o.get_val() <= cutoff)
					break;
			}

			if (s >= len) { // No elements found < cutoff. Throw away sparse
				// list
				Link n = l.next;
				remove(l);
				l = n;
				continue;
			}

			// Find last element that is < cutoff
			for (e = len; e > s; e--) {
				Has_Value o = (Has_Value) l.o.get(e - 1);
				if (o.get_val() <= cutoff)
					break;
			}

			l.end = e + l.start;
			l.start = s + l.start;

			// Remove elements from the start of the list.
			for (int i = 0; i < s; i++) {
				delObj((Has_Value) l.o.remove(0));
			}

			l.realloc(); // This will remove any elements from the end of the
			// list

			l = l.next;
		}
		sanity();
	}

	private Has_Value newObj() {
		if (unusedObj.size() > 0)
			return (Has_Value) unusedObj.removeFirst();
		newObj++;
		return (Has_Value) ((Has_Value) type).clone();
	}

	private void delObj(Has_Value o) {
		unusedObj.addFirst(o);
	}

	// moveFwd(int, Iterate) - find the element at 'int'. iterate.o==null if
	// none there.
	// Iterate position updated.
	Iterate moveFwd(int n, Iterate res) {
		return moveFwd(n, res, false);
	}

	// moveFwd(Iterate) - find the next element. iterate.o==null if at end.
	// Iterate position updated.
	Iterate moveFwd(Iterate res) {
		return moveFwd((res == null ? 0 : res.i + 1), res, false);
	}

	Iterate moveFwd(int n, Iterate res, boolean exact) {
		if (res == null) {
			res = new Iterate();
			res.curr = head;
		}

		res.o = null;
		while (res.curr != null) {
			if (exact && n < res.curr.start)
				return res; // Didn't find an element at n

			if (n < res.curr.end) {
				if (n < res.curr.start)
					n = res.curr.start; // Update n to next available item
				res.i = n;
				res.o = (Has_Value) res.curr.o.get(n - res.curr.start);
				return res;
			}
			res.curr = res.curr.next;
		}
		return res;
	}

	Iterate moveRev(int n, Iterate res) {
		return moveRev(n, res, false);
	}

	Iterate moveRev(Iterate res) {
		int end = (tail == null ? 0 : tail.end);
		return moveRev((res == null ? end : res.i - 1), res, false);
	}

	Iterate moveRev(int n, Iterate res, boolean exact) {
		if (res == null) {
			res = new Iterate();
			res.curr = tail;
		}

		res.o = null;
		while (res.curr != null) {
			if (exact && n >= res.curr.end)
				return res; // Didn't find an element at n

			if (n >= res.curr.start) {
				if (n >= res.curr.end)
					n = res.curr.end - 1;
				res.i = n;
				res.o = (Has_Value) res.curr.o.get(n - res.curr.start);
				return res;
			}
			res.curr = res.curr.prev;
		}
		return res;
	}

	// get the object at interate posSrc + 1. Expected to be in same link.
	Has_Value getNext(Iterate res) {
		if (res.curr == null)
			return null;

		int n = res.i + 1;
		if (n < res.curr.start || n >= res.curr.end)
			return null;
		return (Has_Value) res.curr.o.get(n - res.curr.start);
	}

	// Get the object at interate posSrc - 1. Expected to be in same link.
	Has_Value getPrev(Iterate res) {
		if (res.curr == null)
			return null;

		int n = res.i - 1;
		if (n < res.curr.start || n >= res.curr.end)
			return null;
		return (Has_Value) res.curr.o.get(n - res.curr.start);
	}

	void incAge() {
		for (Link l = head; l != null; l = l.next) {
			l.age++;
		}
	}

	public String toString() {
		StringBuffer r = new StringBuffer();
		for (Link l = head; l != null; l = l.next) {
			r.append(l);
		}
		return r.toString();
	}

	public void sanity() {
		Link p = null;
		for (Link l = head; l != null; l = l.next) {
			Misc.my_assert(l.start >= 0 && l.start < l.end,
					"Insane. start/end index stuffed.");

			if (l.prev != null)
				Misc.my_assert(l == l.prev.next, "Insane 1");
			if (l.next != null)
				Misc.my_assert(l == l.next.prev, "Insane 2");

			p = l;
		}

		if (head != null)
			Misc.my_assert(head.prev == null, "Insane 4");

		if (tail != null)
			Misc.my_assert(tail.next == null, "Insane 5");

		if (p != null)
			Misc.my_assert(p == tail, "Insane 6");
	}

	public void sanity2() {
		Link p = null;
		for (Link l = head; l != null; l = l.next) {
			Misc.my_assert(l.start >= 0 && l.start < l.end,
					"Insane. start/end index stuffed.");
			if (p != null) {
				if (p.end > l.start) {
					System.err.println("Insane.  Regions overlap");
					System.err.println("Sparse = " + this);
					Misc.my_assert(false, "");
				}
			}
			p = l;
		}
	}

	public void display_stats() {
		display_stats("");
	}

	public void display_stats(String pre) {
		long waste = 0;
		long num = 0;
		long numActive = 0;
		for (Link l = head; l != null; l = l.next) {
			numActive += (l.end - l.start);
			waste += l.o.size() - (l.end - l.start);
			num++;
		}
		System.out.println(pre + "Number of active lists = " + num);
		System.out.println(pre + "Total active cells     = " + numActive);
		System.out.println(pre + "Wasted array elements  = " + waste);
		System.out.println(pre + "Total links created    = " + links_created);
		System.out.println(pre + "Total Obj created      = " + newObj);
		// System.out.println(this);
	}
}