<?php
/**
 * Plugin Name:       DeGenPrime
 * Plugin URI:        https://github.com/raw-lab/degenprime
 * Description:       Degenerate PCR primer design from DNA sequences, aligned or not (MAFFT is built in). The Rust engine runs in the visitor's browser as WebAssembly, so sequences are never uploaded. Add [degenprime] to any page.
 * Version:           2.3.0
 * Requires at least: 5.8
 * Requires PHP:      7.4
 * Author:            RAW Lab, UNC Charlotte
 * Author URI:        https://www.rawlab.org
 * License:           Non-commercial; see LICENSE
 * Text Domain:       degenprime
 */

if ( ! defined( 'ABSPATH' ) ) {
	exit;
}

define( 'DEGENPRIME_VERSION', '2.3.0' );

/** Views the app can open on. */
function degenprime_views() {
	return array( 'design', 'results', 'check', 'manual', 'about' );
}

function degenprime_assets_url() {
	return plugins_url( 'assets/', __FILE__ );
}

function degenprime_register_assets() {
	wp_register_script( 'degenprime-app', degenprime_assets_url() . 'degenprime.js', array(), DEGENPRIME_VERSION, true );
	wp_register_style( 'degenprime-embed', plugins_url( 'degenprime-embed.css', __FILE__ ), array(), DEGENPRIME_VERSION );
}
add_action( 'wp_enqueue_scripts', 'degenprime_register_assets' );

/**
 * [degenprime view="design" layout="wide"]
 *
 * view   design | results | check | manual | about   Page the app opens on.
 * layout contained | wide | full                     How far the app may extend past the theme's content column.
 */
function degenprime_shortcode( $atts ) {
	$atts   = shortcode_atts( array( 'view' => 'design', 'layout' => 'wide' ), $atts, 'degenprime' );
	$view   = in_array( $atts['view'], degenprime_views(), true ) ? $atts['view'] : 'design';
	$layout = in_array( $atts['layout'], array( 'contained', 'wide', 'full' ), true ) ? $atts['layout'] : 'wide';
	if ( ! wp_script_is( 'degenprime-app', 'registered' ) ) {
		degenprime_register_assets();
	}
	wp_enqueue_script( 'degenprime-app' );
	wp_enqueue_style( 'degenprime-embed' );
	return sprintf(
		'<div class="degenprime-embed degenprime-layout-%1$s"><degenprime-app data-assets="%2$s" data-view="%3$s" data-host="wordpress"></degenprime-app><noscript>%4$s</noscript></div>',
		esc_attr( $layout ),
		esc_url( degenprime_assets_url() ),
		esc_attr( $view ),
		esc_html__( 'DeGenPrime needs JavaScript and WebAssembly to run.', 'degenprime' )
	);
}
add_shortcode( 'degenprime', 'degenprime_shortcode' );

/** Tools > DeGenPrime: usage notes and a one-click "create page" button. */
function degenprime_admin_menu() {
	add_management_page( 'DeGenPrime', 'DeGenPrime', 'edit_pages', 'degenprime', 'degenprime_admin_page' );
}
add_action( 'admin_menu', 'degenprime_admin_menu' );

function degenprime_admin_page() {
	if ( ! current_user_can( 'edit_pages' ) ) {
		return;
	}
	$created = isset( $_GET['degenprime_created'] ) ? absint( $_GET['degenprime_created'] ) : 0; // phpcs:ignore WordPress.Security.NonceVerification
	?>
	<div class="wrap">
		<h1>DeGenPrime</h1>
		<?php if ( $created ) : ?>
			<div class="notice notice-success"><p>
				<?php esc_html_e( 'Page created.', 'degenprime' ); ?>
				<a href="<?php echo esc_url( get_permalink( $created ) ); ?>"><?php esc_html_e( 'View it', 'degenprime' ); ?></a> |
				<a href="<?php echo esc_url( get_edit_post_link( $created ) ); ?>"><?php esc_html_e( 'Edit it', 'degenprime' ); ?></a>
			</p></div>
		<?php endif; ?>
		<p><?php esc_html_e( 'Add the primer designer to any page or post with this shortcode:', 'degenprime' ); ?> <code>[degenprime]</code></p>
		<table class="widefat striped" style="max-width:760px">
			<thead><tr><th><?php esc_html_e( 'Attribute', 'degenprime' ); ?></th><th><?php esc_html_e( 'Values', 'degenprime' ); ?></th><th><?php esc_html_e( 'Default', 'degenprime' ); ?></th></tr></thead>
			<tbody>
				<tr><td><code>view</code></td><td>design, results, check, manual, about</td><td>design</td></tr>
				<tr><td><code>layout</code></td><td>contained, wide (up to 1240 px), full (edge to edge)</td><td>wide</td></tr>
			</tbody>
		</table>
		<h2><?php esc_html_e( 'Create a page in one click', 'degenprime' ); ?></h2>
		<form method="post" action="<?php echo esc_url( admin_url( 'admin-post.php' ) ); ?>">
			<input type="hidden" name="action" value="degenprime_create_page">
			<?php wp_nonce_field( 'degenprime_create_page' ); ?>
			<?php submit_button( __( 'Create the "Primer design" page', 'degenprime' ), 'primary', 'submit', false ); ?>
		</form>
		<h2><?php esc_html_e( 'If the app reports that the engine could not start', 'degenprime' ); ?></h2>
		<ul style="list-style:disc;padding-left:20px;max-width:760px">
			<li><?php esc_html_e( 'A Content-Security-Policy header (often set by a security plugin) must allow \'wasm-unsafe-eval\' in script-src and blob: in worker-src.', 'degenprime' ); ?></li>
			<li><?php esc_html_e( 'Script optimisation plugins should not combine or delay degenprime.js.', 'degenprime' ); ?></li>
			<li><?php esc_html_e( 'The bundled assets/.htaccess serves .wasm files with the right type on Apache; on nginx add "application/wasm wasm;" to mime.types (the app also works without it, just a little slower to start).', 'degenprime' ); ?></li>
		</ul>
	</div>
	<?php
}

function degenprime_create_page() {
	if ( ! current_user_can( 'publish_pages' ) ) {
		wp_die( esc_html__( 'You are not allowed to publish pages.', 'degenprime' ) );
	}
	check_admin_referer( 'degenprime_create_page' );
	$page_id = wp_insert_post(
		array(
			'post_title'   => __( 'Primer design', 'degenprime' ),
			'post_content' => "<!-- wp:shortcode -->\n[degenprime layout=\"wide\"]\n<!-- /wp:shortcode -->",
			'post_status'  => 'publish',
			'post_type'    => 'page',
		),
		true
	);
	if ( is_wp_error( $page_id ) ) {
		wp_die( esc_html( $page_id->get_error_message() ) );
	}
	wp_safe_redirect( add_query_arg( 'degenprime_created', $page_id, admin_url( 'tools.php?page=degenprime' ) ) );
	exit;
}
add_action( 'admin_post_degenprime_create_page', 'degenprime_create_page' );

function degenprime_action_links( $links ) {
	array_unshift( $links, '<a href="' . esc_url( admin_url( 'tools.php?page=degenprime' ) ) . '">' . esc_html__( 'Set up', 'degenprime' ) . '</a>' );
	return $links;
}
add_filter( 'plugin_action_links_' . plugin_basename( __FILE__ ), 'degenprime_action_links' );
