package japsadev.bio.phylo;

import java.text.Normalizer;
import java.text.Normalizer.Form;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Slug{
	private static final Pattern NONLATIN = Pattern.compile("[^\\w-]");
	private static final Pattern WHITESPACE = Pattern.compile("[\\s_]");
	public static String toSlug(String input) {
	  Matcher matcher = WHITESPACE.matcher(input);
	 String nowhitespace = WHITESPACE.matcher(input).replaceAll("_");
	  String normalized = Normalizer.normalize(nowhitespace, Form.NFD);
	  String slug = NONLATIN.matcher(normalized).replaceAll("");
	  return slug.toLowerCase(Locale.ENGLISH);
	}
}