package japsafuture.stream;

import static java.lang.Thread.currentThread;
import static java.util.concurrent.Executors.newSingleThreadExecutor;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Flow.Publisher;
import java.util.concurrent.Flow.Subscriber;
import java.util.concurrent.Flow.Subscription;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

public class MyPublisher implements Publisher<Integer> {

    private static final String LOG_MESSAGE_FORMAT = "Publisher >> [%s] %s%n";

    final ExecutorService executor = Executors.newFixedThreadPool(4);
    private List<MySubscription> subscriptions = Collections.synchronizedList(new ArrayList<MySubscription>());

    private final CompletableFuture<Void> terminated = new CompletableFuture<>();

    @Override
    public void subscribe(Subscriber<? super Integer> subscriber) {
        MySubscription subscription = new MySubscription(subscriber, executor);

        subscriptions.add(subscription);

        subscriber.onSubscribe(subscription);
    }

    public void waitUntilTerminated() throws InterruptedException {
        try {
            terminated.get();
        } catch (ExecutionException e) {
            System.out.println(e);
        }
    }

    private class MySubscription implements Subscription {

        private final ExecutorService executor;

        private Subscriber<? super Integer> subscriber;
        private final AtomicInteger value;
        private AtomicBoolean isCanceled;

        public MySubscription(Subscriber<? super Integer> subscriber, ExecutorService executor) {
            this.subscriber = subscriber;
            this.executor = executor;

            value = new AtomicInteger();
            isCanceled = new AtomicBoolean(false);
        }

        @Override
        public void request(long n) {
            if (isCanceled.get())
                return;

            if (n < 0)
                executor.execute(() -> subscriber.onError(new IllegalArgumentException()));
            else
                publishItems(n);
        }

        @Override
        public void cancel() {
            isCanceled.set(true);

            synchronized (subscriptions) {
                subscriptions.remove(this);
                if (subscriptions.size() == 0)
                    shutdown();
            }
        }

        private void publishItems(long n) {
            for (int i = 0; i < n; i++) {

                executor.execute(() -> {
                    int v = value.incrementAndGet();
                    log("publish item: [" + v + "] ...");
                    subscriber.onNext(v);
                });
            }
        }

        private void shutdown() {
            log("Shut down executor...");
            executor.shutdown();
            newSingleThreadExecutor().submit(() -> {

                log("Shutdown complete.");
                terminated.complete(null);
            });
        }

    }

    private void log(String message, Object... args) {
        String fullMessage = String.format(LOG_MESSAGE_FORMAT, currentThread().getName(), message);

        System.out.printf(fullMessage, args);
    }
}